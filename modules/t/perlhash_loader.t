use strict;
use warnings;
use Test::More;
use Test::Differences;
use Test::Exception;
use Bio::EnsEMBL::Test::TestUtils qw/warns_like/;

require_ok('Bio::EnsEMBL::Registry');
require_ok('Bio::EnsEMBL::Registry::Loader::PerlHash');

my $registry = 'Bio::EnsEMBL::Registry';
my $verbose = 0;
my $loader = Bio::EnsEMBL::Registry::Loader::PerlHash->new($registry, $verbose);

#### Check for empty attempts at loading a registry

eq_or_diff($loader->serialise_registry(), {adaptors => [], aliases => {}}, 'An empty registry means a minimal perl data structure');
$loader->load_registry({aliases => {}});
eq_or_diff($loader->serialise_registry(), {adaptors => [], aliases => {}}, 'No adaptors in the submitted hash is tolerated');
$loader->load_registry({adaptors => []});
eq_or_diff($loader->serialise_registry(), {adaptors => [], aliases => {}}, 'No aliases in the submitted hash is tolerated');

#### Check that we die if we give bad data structures

dies_ok { $loader->load_registry({aliases => []})} 'Making aliases an array throws an exception';
dies_ok { $loader->load_registry({adaptors => {}})} 'Making adaptors a hash throws an exception';

#### Check actually loading some real data that makes sense

my $perl_hash = {
  adaptors => [
    {driver => 'mysql', dbname => 'fake', host => 'localhost', port => 3306, user => 'user', pass => 'pass', 
      species => 'human', group => 'core'},
    {driver => 'mysql', dbname => 'fake', host => 'localhost', port => 3306, user => 'user', pass => 'pass', 
      species => 'mouse', group => 'core'},
    {driver => 'mysql', dbname => 'fake', host => 'localhost', port => 3306, user => 'user', pass => 'pass', 
      species => 'ecoli', group => 'core', multispecies_db => 1, species_id => 20},
  ],
  aliases => {
    human => [qw/9606 homer/],
    mouse => [qw/mice squeek/],
  }
};

$loader->load_registry($perl_hash);
my $produced_perl_hash = $loader->serialise_registry();

is(scalar(@{$registry->get_all_DBAdaptors()}), 3, 'Expected number of DBAdaptors in the Registry');
eq_or_diff($produced_perl_hash, $perl_hash, 'Checking the data roundtrips from PerlHash');

$registry->clear();

#### Check for bad groups

warns_like {
  $loader->load_registry({adaptors => [ {group => 'wibble', species => 'bad'}, {group => 'wibble', species => 'badagain'}]});
} qr/wibble .+ known .+ module/xms, 'Checking that an unknown group means an adequate warning';
$registry->clear();

#### Check that we RESPECT the no cache flag

my $no_cache = 1;
my $no_cache_loader = Bio::EnsEMBL::Registry::Loader::PerlHash->new($registry, 0, $no_cache);
$no_cache_loader->load_registry($perl_hash);
ok($registry->get_DBAdaptor('human', 'core')->no_cache(), 'Check no_cache was propagated through');
$registry->clear();

done_testing();