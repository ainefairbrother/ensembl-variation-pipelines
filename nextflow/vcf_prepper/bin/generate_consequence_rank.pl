#!/usr/bin/env perl

#################################################################
#
# Prepare a json config with variation consequnce and their rank
# Depend on the rank configured in ensembl-variation
#
#################################################################

use strict;
use warnings;
use JSON;
use Getopt::Long;

use Bio::EnsEMBL::Variation::Utils::Constants;

# parse cli parameters
my $config = {};
GetOptions(
  $config,
  
  'help',                     # displays help message
  'output_file|o=s',          # output file name
) or die "ERROR: Failed to parse command-line flags\n";
&usage && exit(0) if $config->{help};

# help message
sub usage {
  my $usage =<<END;
Genrate a json file with variation consequence and their ranks using -
Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES

Usage:
  perl generate_consequence_rank.pl [--output_file]

  Basic options
  =============

  --help                 Display this message and quit
  -o | --output_file     Output file, by default "variation_consequnce_rank.json"
END

    print $usage;
}

# getting the overlap consequece object from ensembl-variation repo
my %all_cons = %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;

# parse all_cons hash object and create a hash with only rank values
my $cons_rank = {};
foreach my $cons (keys %all_cons){
  $cons_rank->{$cons} = $all_cons{$cons}->rank;
}

# convert hash to pretty json
my $cons_rank_json = JSON->new->ascii->pretty->encode($cons_rank);

# dump the json to a output file
my $output_file = $config->{"output_file"} || "variation_consequnce_rank.json";
open(my $output_fh, ">", $output_file) or die "ERROR: Failed to open $output_file - $!\n";
print "Printing variation consequnce ranks ...\n";
print $output_fh $cons_rank_json;
close($output_fh);