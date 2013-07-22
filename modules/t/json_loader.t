use strict;
use warnings;
use Test::More;
use Test::Differences;
use Test::Exception;

require_ok('Bio::EnsEMBL::Registry');
require_ok('Bio::EnsEMBL::Registry::Loader::JSON');

my $registry = 'Bio::EnsEMBL::Registry';
my $verbose = 0;
my $loader = Bio::EnsEMBL::Registry::Loader::JSON->new($registry, $verbose);

#### Check for empty attempts at loading a registry

# JSON uses pretty print 3 spacing for some reason
my $standard_empty_json = <<JSON;
{
   "adaptors" : [],
   "aliases" : {}
}
JSON
eq_or_diff($loader->serialise_registry(), $standard_empty_json, 'An empty registry means a minimal perl data structure');
$loader->load_registry('{"aliases" : {}}');
is($loader->serialise_registry(), $standard_empty_json, 'No adaptors in the submitted JSON is tolerated');
$loader->load_registry('{"adaptors" : []}');
is($loader->serialise_registry(), $standard_empty_json, 'No aliases in the submitted JSON is tolerated');

#### Check that we die if we give bad data structures

dies_ok { $loader->load_registry(undef) } 'Initalising with an undef means death';
dies_ok { $loader->load_registry(q{}) } 'Initalising without a JSON document means death';

#### Check actually loading some real data that makes sense

# Check basic JSON
my $json = <<JSON;
{
   "adaptors" : [
      {
         "dbname" : "fake",
         "driver" : "mysql",
         "group" : "core",
         "host" : "localhost",
         "pass" : "pass",
         "port" : 3306,
         "species" : "human",
         "user" : "user"
      },
      {
         "dbname" : "fake",
         "driver" : "mysql",
         "group" : "core",
         "host" : "localhost",
         "pass" : "pass",
         "port" : 3306,
         "species" : "mouse",
         "user" : "user"
      },
      {
         "dbname" : "fake",
         "driver" : "mysql",
         "group" : "core",
         "host" : "localhost",
         "multispecies_db" : 1,
         "pass" : "pass",
         "port" : 3306,
         "species" : "ecoli",
         "species_id" : 20,
         "user" : "user"
      }
   ],
   "aliases" : {
      "human" : [
         "9606",
         "homer"
      ],
      "mouse" : [
         "mice"
      ]
   }
}
JSON

$loader->load_registry($json);
my $produced_json = $loader->serialise_registry();

is(scalar(@{$registry->get_all_DBAdaptors()}), 3, 'Expected number of DBAdaptors in the Registry');
eq_or_diff($produced_json, $json, 'Checking the data roundtrips from JSON');

$registry->clear();

## Check if we try to insert some inline comments and the alike
my $human_json = $json;
$human_json =~ s/"homer"/"homer", #this is no-longer strict JSON/;
$loader->load_registry($human_json);
like($human_json, qr/, #this is/, 'Checking the JSON was modified as expected');
eq_or_diff($loader->serialise_registry(), $json, 'Checking the data roundtrips from more user generated JSON');

done_testing();