use strict;
use warnings;
use Test::More;
use Test::Differences;
use Test::Exception;
use File::Temp qw/tempfile/;
use Bio::EnsEMBL::Test::TestUtils qw/warns_like/;

require_ok('Bio::EnsEMBL::Registry');
require_ok('Bio::EnsEMBL::Registry::Loader::ConfigIni');

my $registry = 'Bio::EnsEMBL::Registry';
my $verbose = 0;
my $loader = Bio::EnsEMBL::Registry::Loader::ConfigIni->new($registry, $verbose);

sub write_inifile {
  my ($input) = @_;
  my ($fh, $filename) = tempfile();
  print $fh $input;
  close $fh;
  return ($fh, $filename); #keep $fh alive otherwise the file goes
}

#### Check for empty attempts at loading a registry
{
  my ($fh, $filename) = write_inifile('');
  $loader->load_registry($filename);
  is($loader->serialise_registry()->Sections(), 0, 'No adaptors in the submitted INI is tolerated');
}
{
  warns_like { $loader->load_registry('/a/random/location')} qr/No such file or directory/, 'Bad location means no config';
  is($loader->serialise_registry()->Sections(), 0, 'No adaptors in the submitted INI is tolerated');
}
{
  my ($fh, $filename) = write_inifile('[default]');
  $loader->load_registry($filename);
  is($loader->serialise_registry()->Sections(), 0, 'Just a [default] section does nothing');
}

#### Check actually loading some real data that makes sense

# Check basic INI file
my $ini = "[default]
host=somewhere

[human_core]
host=localhost
port=3306
user=user
pass=pass
dbname=db
species=human
group=core
alias=<<ALIAS
9606\r\nhomer\ntest
ALIAS

[ecoli_core]
port=3306
user=user
pass=pass
dbname=db
species=ecoli
group=core
multispecies_db=1
species_id=20
";

{
  my ($fh, $filename) = write_inifile($ini);
  $loader->load_registry($filename);
  my $produced_ini = $loader->serialise_registry();
  is($registry->get_DBAdaptor_count(), 2, 'Expected number of DBAdaptors in the Registry');
  is($produced_ini->Sections(), 2, 'Checking available section count');
  
  #Check human
  eq_or_diff([$produced_ini->val('human_core', 'alias')], ['test','9606','homer'], 'Check aliases are OK and can parse windows line terminators');
  
  #Check ecoli
  is($produced_ini->val('ecoli_core', 'host'), 'somewhere', 'Check default host came through');
  is($produced_ini->val('ecoli_core', 'multispecies_db'), 1, 'Check multispecies_db came through');
  is($produced_ini->val('ecoli_core', 'species_id'), 20, 'Check species id came through');
}
$registry->clear();

done_testing();