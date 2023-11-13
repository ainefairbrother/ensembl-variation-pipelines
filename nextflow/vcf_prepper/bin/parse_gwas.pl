#!/usr/bin/env perl

#################################################################
#
# Parse GWAS tsv file into a VCF
# (Not directly used by the pipeline currently) 
#
#################################################################

use Path::Tiny qw(path);
use Storable qw(dclone);

use Bio::EnsEMBL::Registry;

my ($input_file, $workdir) = @ARGV;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($workdir . '/ensembl.registry');
$registry->add_alias('homo_sapiens', 'human');

my $va = $registry->get_adaptor('human', 'variation', 'variation');

sub get_vfs_from_id {
  my ($id) = @_;
  
  return [] unless defined $va;
  
  my $v = $va->fetch_by_name($id);
  return [] unless defined $v;
  
  my $locations = [];
  foreach my $vf (@{ $v->get_all_VariationFeatures() }) {
    
    my $seq = $vf->seq_region_name();
    my $start = $vf->seq_region_start();
    my $end = $vf->seq_region_end();
    my $ref = $vf->ref_allele_string();
    
    if($vf->ref_allele_string =~ /-/) {
      # convert to vcf format to compare the alt alleles
      my $convert_to_vcf = $vf->to_VCF_record;
      $start = ${$convert_to_vcf}[1];
      $ref = ${$convert_to_vcf}[3];
    }
    
    next unless $ref =~ /[ACGTN]+/;
    
    my $location = {
      "seq"   => $seq,
      "start" => $start,
      "end"   => $end,
      "ref"   => $ref
    };
    
    push @{ $locations }, $location;
  }
  
  return $locations;
}

sub parse_input_file {
  my ($input_file) = @_;

  my $input_filename = path($input_file)->basename();
  my (%headers, $variants);

  my ($err_file, $err_FH);
  $err_file = $workdir . '/log_' . $input_filename . ".err";
  open($err_FH, ">", $err_file) || die ("Could not open $err_file for writing: $!\n");

  # Open the input file for reading
  if($input_file =~ /gz$/) {
    open(IN, "zcat " . $input_file . " |") || die ("Could not open $input_file for reading: $!\n");
  }
  else {
    open(IN, '<', $input_file) || die ("Could not open $input_file for reading: $!\n");
  }

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    my @row_data = split(/\t/,$_);

    # header
    if(/^DATE\s+ADDED\s+TO\s+CATALOG/) {
      $headers{uc($row_data[$_])} = $_ for 0..$#row_data;
    }
    else {
      die ("ERROR: Could not find header data\n") unless %headers;

      my %content;
      $content{$_} = $row_data[$headers{$_}] for keys %headers;

      my $rs_risk_allele = ($content{'STRONGEST SNP-RISK ALLELE'} =~ /\?/) ? '' : $content{'STRONGEST SNP-RISK ALLELE'};
      my $rs_id          = $content{'SNPS'};
      

      # Parse the ids
      my @ids;
      $rs_id ||= "";
      while ($rs_id =~ m/(rs[0-9]+)/g) {
        push(@ids, $1);
      }

      print $err_FH "WARNING: Could not parse any rsIds from string '$rs_id'\n" if (!scalar(@ids));
      next if (!scalar(@ids));

      map {
        my $id = $_;
        my $t_data = {};
        
        my ($risk_allele, $risk_allele_with_id);
        map {
          if ($_ =~ /$id/) {
            $risk_allele_with_id = $_;
            $risk_allele = ( split("-", $risk_allele_with_id) )[1];
          }
        } split(/[;,x]/, $rs_risk_allele);
        
        unless (defined $variants->{$risk_allele_with_id}){
          if ($risk_allele =~ /[ATCGN]+/){
            my $vfs = get_vfs_from_id($id);
            $variants->{$risk_allele_with_id} = {
              "id" => $id,
              "allele" => $risk_allele,
              "vfs" => $vfs
            }
          }
        }
      } @ids;
    }
  }
  close(IN);
  close($err_FH);

  return $variants;
}

sub create_input_vcf {
  my ($data, $input_file) = @_;
  
  open (my $input_FH, ">", $input_file) || die ("Could not open $input_file for writing: $!\n");
  
  foreach (keys %{ $data }){
    my $var = $data->{$_};
    my $vf = $var->{"vfs"};
    
    foreach (@{ $vf }){
      my $line = join("\t", (
          $_->{"seq"}, $_->{"start"}, $var->{"id"}, $_->{"ref"}, $var->{"allele"}, ".", "."
      ));
        
      print $input_FH $line . "\n";
    }
  }
  
  close($input_FH);
}

my $result = parse_input_file($input_file);
my $vcf_file = $workdir . "/input_" . path($input_file)->basename(".tsv") . ".vcf";
create_input_vcf($result, $vcf_file);
system("sort -k1,1 -k2,2n $vcf_file") == 0 
  || die ("Failed to sort $vcf_file:\n$?\n");
