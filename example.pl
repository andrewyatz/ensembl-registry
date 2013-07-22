use Bio::EnsEMBL::Registry;

#Bio::EnsEMBL::Registry->load_registry_from_db(-HOST => 'ensembldb.ensembl.org', -USER => 'anonymous', -DB_VERSION => $version, -VERBOSE => 1);
#my $dba =Bio::EnsEMBL::Registry->get_DBAdaptor('homo_sapiens', 'core');
#warn $_ for DBI->data_sources('mysql', {host => 'ensembldb.ensembl.org', port => 5306, user => 'anonymous'}); 
#warn '';
#warn '';
#warn '';
#warn $_ for @{$dba->dbc->sql_helper->execute_simple(-SQL => 'SHOW DATABASES LIKE ?', -PARAMS => ["%$version%"])}; 

$| =1;
#use re 'debug';
use strict;
use warnings;
use Config::IniFiles;
use Bio::EnsEMBL::Registry::Storage;
use Bio::EnsEMBL::Registry::Loader::ConfigIni;
use Bio::EnsEMBL::Registry::Loader::JSON;
use Bio::EnsEMBL::Registry::Loader::Database;
use Bio::EnsEMBL::Utils::IO qw/slurp/;
use IO::String;

my $version = 63;

my $verbose = 0;

use Benchmark qw/cmpthese/;

cmpthese(100, {
#  old => sub {
#    Bio::EnsEMBL::Registry->load_registry_from_db(-HOST => 'ensembldb.ensembl.org', -USER => 'anonymous', -DB_VERSION => $version, -VERBOSE => $verbose);
#    Bio::EnsEMBL::Registry->clear();
#  },
#  new => sub {
#    my $storage = Bio::EnsEMBL::Registry::Storage->new();
#    my $db = Bio::EnsEMBL::Registry::Loader::Database->new(-REGISTRY_STORAGE => $storage);
#    $db->load_from_db(-HOST => 'ensembldb.ensembl.org', -USER => 'anonymous', -DB_VERSION => $version, -VERBOSE => $verbose);
#  },
  config_ini => sub {
    my $storage = Bio::EnsEMBL::Registry::Storage->new();
    my $cfg_ini = Bio::EnsEMBL::Registry::Loader::ConfigIni->new(-REGISTRY_STORAGE => $storage);
    my $cfg = Config::IniFiles->new(-file => 'reg.ini');
    $cfg_ini->load_from(-CFG => $cfg);
  },
  json => sub {
    my $storage = Bio::EnsEMBL::Registry::Storage->new();
    my $json = Bio::EnsEMBL::Registry::Loader::JSON->new(-REGISTRY_STORAGE => $storage);
    my $json_scalar = slurp('reg.json');
    $json->load_from(-JSON => $json_scalar);
  }
});

#my $storage = Bio::EnsEMBL::Registry::Storage->new();
#my $db = Bio::EnsEMBL::Registry::Loader::Database->new(-REGISTRY_STORAGE => $storage);
#my $json = Bio::EnsEMBL::Registry::Loader::JSON->new(-REGISTRY_STORAGE => $storage);
#$db->load_from_db(-HOST => 'ensembldb.ensembl.org', -USER => 'anonymous', -DB_VERSION => $version, -VERBOSE => $verbose);
#
#my $json_string = $json->convert_to();
#warn $json_string;

#my $cfg_ini = Bio::EnsEMBL::Registry::Loader::ConfigIni->new();
#$db->load_from_db(-REGISTRY_STORAGE => $storage, -HOST => 'mysql.ebi.ac.uk', -PORT => 4157, -USER => 'anonymous', -DB_VERSION => $version, -VERBOSE => $verbose);



#my $cfg = $cfg_ini->convert_to($storage);
#$cfg->WriteConfig('reg.ini');

#my $ref = q{};
#my $io = IO::String->new(\$ref);
#select($io);
#my $cfg = $cfg_ini->convert_to_ini($storage);
#select();
#print STDOUT $ref, "\n";
#$cfg->OutputConfig();
#warn 'ello';

#my $db_name = 'ailuropoda_melanoleuca_core_63_1';
#my $re = $db->_patterns()->{FULL}->{core}->{single};
#
#warn $re;
#
#if($db_name =~ $re) {
#  print "original one\n";
#}

#my $other_reg = qr/^ ([a-z]+ _ [a-z0-9]+ ) _core (?:_\d+)? _ (\d+) _\w+$/xms;
#my $new_reg = qr/(?msx-i:^ ((?msx-i: [a-z]+ _ [a-z0-9]+ ))   (?msx-i:_core           ) (?msx-i:(?-xism: (?:_\d+)? ) _ ((?-xism: \d+ )) (?msx-i: _\w+ )) $)/;
#my $new_reg = qr/(?msx-i:^ ((?msx-i: [a-z]+ _ [a-z0-9]+ ))   (?msx-i:_core           ) (?msx-i: (?:_\d+)?  _ ((?-xism: \d+ )) (?msx-i: _\w+ )) $)/;
#my $new_reg = qr/(?msx-i:^ ((?msx-i: [a-z]+ _ [a-z0-9]+ ))   (?msx-i:_core           ) (?msx-i:(?:_\d+)? _ (\d+) _ \w+) $)/;
#
#if($db_name =~ $new_reg) {
#  print "hello from the new one\n";
#};
#
#if($db_name =~ $other_reg) {
#  print "hello from the other one\n";
#};
#
#print "Other: $other_reg \n";
#print "New: $new_reg \n";

#$db->load_from_db(
#  -REGISTRY_STORAGE => $storage,
#  -HOST => 'ensembldb.ensembl.org', -USER => 'anonymous', -DB_VERSION => $version, -VERBOSE => 1
#);