use strict;
use warnings;
use Test::More;
use Test::Differences;
use Test::Exception;
use File::Temp qw/tempfile/;

require_ok('Bio::EnsEMBL::Registry');
require_ok('Bio::EnsEMBL::Registry::Loader::Perl');

my $registry = 'Bio::EnsEMBL::Registry';
my $verbose = 0;
my $loader = Bio::EnsEMBL::Registry::Loader::Perl->new($registry, $verbose);

sub write_tmp {
  my ($input) = @_;
  my ($fh, $filename) = tempfile();
  print $fh $input;
  close $fh;
  return ($fh, $filename); #keep $fh alive otherwise the file goes
}

{
  dies_ok { $loader->load_registry('/some/random/path') } 'A non-locatable file should fail';
}

{
  my ($fh, $filename) = write_tmp('');
  dies_ok { $loader->load_registry($filename) } 'A non-true returning eval should fail';
}

{
  my ($fh, $filename) = write_tmp('1;');
  $loader->load_registry($filename);
  is($registry->get_DBAdaptor_count(), 0, 'Expected number of DBAdaptors in the Registry');
}

{
  my ($fh, $filename) = write_tmp('1;');
  $loader->load_registry($filename);
  is($registry->get_DBAdaptor_count(), 0, 'Expected number of DBAdaptors in the Registry');
}

{
  my ($fh, $filename) = write_tmp(<<REG);
use Bio::EnsEMBL::DBSQL::DBAdaptor
Bio::EnsEMBL::DBSQL::DBAdaptor->new(-HOST => 'host', -USER => 'user', -SPECIES => 'human', -GROUP => 'core', -dbname => 'db');
1;
REG
  $loader->load_registry($filename);
  is($registry->get_DBAdaptor_count(), 1, 'Expected number of DBAdaptors in the Registry');
  is($registry->get_DBAdaptor('human','core')->dbc->host(), 'host', 'Expected host in human DBAdaptor');
}

done_testing();