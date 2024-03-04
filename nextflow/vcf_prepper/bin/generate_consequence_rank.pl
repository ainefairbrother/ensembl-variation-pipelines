#!/usr/bin/env perl

=head1 LICENSE
 
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
 
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
 
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
 
=cut

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