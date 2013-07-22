package Bio::EnsEMBL::Registry::Loader::Perl;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Registry::Loader::Base/;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::IO qw/slurp/;

sub load_registry {
  my ($self, $location) = @_;
  my $contents = slurp($location);
  my $test_eval = eval $contents; ## no critic
  eval { require($location) };
  throw "Could not bring in Perl configuration from '$location': $@" unless $test_eval;
  return;
}

1;