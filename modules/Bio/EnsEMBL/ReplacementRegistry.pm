=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Registry

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  my $registry = 'Bio::EnsEMBL::Registry';

  $registry->load_all("configuration_file");

  $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

=head1 DESCRIPTION

All Adaptors are stored/registered using this module. This module should
then be used to get the adaptors needed.

The registry can be loaded from a configuration file using the load_all
method.

If a filename is passed to load_all then this is used.  Else if the
environment variable ENSEMBL_REGISTRY is set to the name on an existing
configuration file, then this is used.  Else if the file .ensembl_init
in your home directory exist, it is used.

For the Web server ENSEMBL_REGISTRY should be set in SiteDefs.pm.  This
will then be passed on to load_all.


The registry can also be loaded via the method load_registry_from_db
which given a database host will load the latest versions of the Ensembl
databases from it.

The four types of registries are for db adaptors, dba adaptors, dna
adaptors and the standard type.

=head2 db

These are registries for backwards compatibility and enable the
subroutines to add other adaptors to connections.

e.g. get_all_db_adaptors, get_db_adaptor, add_db_adaptor,
remove_db_adaptor are the old DBAdaptor subroutines which are now
redirected to the Registry.

So if before we had

  my $sfa = $self->adaptor()->db()->get_db_adaptor('blast');

We now want to change this to

  my $sfa =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "blast" );


=head2 DBA

These are the stores for the DBAdaptors

The Registry will create all the DBConnections needed now if you set up
the configuration correctly. So instead of the old commands like

  my $db           = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);
  my $exon_adaptor = $db->get_ExonAdaptor;

we should now have just

  my $exon_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "exon" );


=head2 DNA

This is an internal Registry and allows the configuration of a dnadb.
An example here is to set the est database to get its dna data from the
core database.

  ## set the est db to use the core for getting dna data.
  # Bio::EnsEMBL::Utils::ConfigRegistry->dnadb_add( "Homo Sapiens",
  #   "core", "Homo Sapiens", "est" );


=head2 adaptors

This is the registry for all the general types of adaptors like
GeneAdaptor, ExonAdaptor, Slice Adaptor etc.

These are accessed by the get_adaptor subroutine i.e.

  my $exon_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "exon" );

=head1 METHODS

=cut

package Bio::EnsEMBL::ReplacementRegistry;
use strict;
use warnings;

use Scalar::Util qw/blessed/;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::ApiVersion qw/software_version/;

# Bring in the loader objects
use Bio::EnsEMBL::Registry::Loader::Base;
use Bio::EnsEMBL::Registry::Loader::ConfigIni;
use Bio::EnsEMBL::Registry::Loader::Database;
use Bio::EnsEMBL::Registry::Loader::DatabaseUrl;
use Bio::EnsEMBL::Registry::Loader::JSON;
use Bio::EnsEMBL::Registry::Loader::Perl;

#Helper objects
use Bio::EnsEMBL::Utils::StableIdsDBAdaptorLookup;

#Declare the global reigstry. Only access this via $self->storage()
use vars qw(%registry_register);

=head2 storage

  Description : Provides access to the most relevant area of Registry 
                storage. If it is given an unblessed reference then
                the global registry will be returned. Otherwise the
                passed in reference is returned since this is the
                appropriate storage layer
  Returntype  : HashRef to be used for Registry operations
  Example     : my $storage = 'Bio::EnsEMBL::Registry'->storage();  #namespace == global storage
                my $storage = Bio::EnsEMBL::Registry->storage();    #bare class == global storage
                
                my $registry = Bio::EnsEMBL::Registry->new();
                my $storage = $registry->storage();                 #blessed reference == non-global storage

=cut

sub storage {
  my ($self) = @_;
  return blessed($self) ? $self : \%registry_register;
}

=head2 new

  Description : Creates a new Registry. This will be a non-global
                registry.
  Returntype  : Bio::EnsEMBL::Registry blessed instance
  Example     : my $registry = Bio::EnsEMBL::Registry->new();

=cut

sub new {
  my ($class) = @_;
  $class = ref($class) || $class;
  return bless({}, $class);
}

# This is a map from group names to Ensembl DB adaptors.  Used by
# load_all() and reset_DBAdaptor().
my %group2adaptor = (
      'blast'      => 'Bio::EnsEMBL::External::BlastAdaptor',
      'compara'    => 'Bio::EnsEMBL::Compara::DBSQL::DBAdaptor',
      'core'       => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'estgene'    => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'funcgen'    => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
      'regulation' => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
      'haplotype' => 'Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor',
      'hive'      => 'Bio::EnsEMBL::Hive::DBSQL::DBAdaptor',
      'ontology'  => 'Bio::EnsEMBL::DBSQL::OntologyDBAdaptor',
      'otherfeatures' => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'pipeline'      => 'Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor',
      'snp'       => 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor',
      'stable_ids' => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'variation' => 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
      'vega'      => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
      'vega_update' => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
);


=head2 load_all

 Will load the registry with the configuration file which is
 obtained from the first in the following and in that order.

  1) If an argument is passed to this method, this is used as the
     name of the configuration file to read.

  2) If the environment variable ENSEMBL_REGISTRY is set, this is
     used as the name of the configuration file to read.

  3) If the file .ensembl_init exist in the home directory, it is
     used as the configuration file.

  Arg [1]    : (optional) string
               Name of file to load the registry from.

  Arg [2]    : (optional) integer
               If not 0, will print out all information.

  Arg [3]    : (optional) integer
               If not 0, the database connection will not be
               cleared, if 0 or if not set the database connections
               will be cleared (this is the default).

  Arg [4]:     (optional) boolean
               This option will turn off caching for slice features,
               so, every time a set of features is retrieved,
               they will come from the database instead of the
               cache.  This option is only recommended for advanced
               users, specially if you need to store and retrieve
               features.  It might reduce performance when querying
               the database if not used properly.  If in doubt, do
               not use it or ask in the developer mailing list.

  Example    : Bio::EnsEMBL::Registry->load_all();
  Returntype : Int count of the DBAdaptor instances which can be found in the 
               registry due to this method being called. Will never be negative
  Exceptions : none
  Status     : Stable

=cut

sub load_all {
  my ($self, $config_file, $verbose, $no_clear, $no_cache ) = @_;

  if ( !defined($config_file) ) {
    if ( defined( $ENV{ENSEMBL_REGISTRY} ) ) {
      $config_file = $ENV{ENSEMBL_REGISTRY};
    } 
    elsif ( defined( $ENV{HOME} ) ) {
      $config_file = $ENV{HOME} . "/.ensembl_init";
    }
  }

  $verbose  ||= 0;
  $no_clear ||= 0;
  $no_cache ||= 0;
  
  my $storage = $self->storage();
  
  my $original_count = $self->get_DBAdaptor_count();

  if ( !defined($config_file) ) {
    print STDERR "No default registry configuration to load.\n" if $verbose;
  } 
  elsif ( ! -e $config_file ) {
    printf( STDERR "Configuration file '%s' does not exist. Registry configuration not loaded.\n", $config_file ) if $verbose;
  } 
  else {
    if ( defined( $storage->{'seen'} ) ) {
      if ( !$no_clear ) {
        print STDERR "Clearing previously loaded registry configuration\n" if $verbose;
        $self->clear();
      }
    }
    $storage->{'seen'} = 1;

    printf( STDERR "Loading registry configuration from '%s'.\n", $config_file) if $verbose;
    
    #If it ends in a .ini then assume ConfigIni
    if($config_file =~ /\.ini$/) {
      my $loader = Bio::EnsEMBL::Registry::Loader::ConfigIni->new($self, $verbose, $no_cache);
      $loader->load_registry($config_file);
    }
    elsif($config_file =~ /\.json$/) {
      my $loader = Bio::EnsEMBL::Registry::Loader::JSON->new($self, $verbose, $no_cache);
      $loader->load_registry($config_file);
    }
    #Otherwise assume it is a Perl file to be evaluated
    else {
      my $perl_eval_loader = Bio::EnsEMBL::Registry::Loader::Perl->new($self, $verbose, $no_cache);
      $perl_eval_loader->load_registry($config_file);
    }
  }
  
  my $count = $self->get_DBAdaptor_count() - $original_count;
  return $count >= 0 ? $count : 0; 
} ## end sub load_all

=head2 clear

 Will clear the registry and disconnect from all databases.

  Example    : Bio::EnsEMBL::Registry->clear();
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub clear {
  my ($self) = @_;
  my $storage = $self->storage();
  foreach my $dba (@{$storage->{'_DBA'}}){
    if($dba->dbc->connected){
      $dba->dbc->db_handle->disconnect();
    }
  }
  %{$storage} = ();
  return;
}

#
# db adaptors. (for backwards compatibility)
#

=head2 add_db

  Arg [1]    : db (DBAdaptor) to add adaptor to.
  Arg [2]    : name of the name to add the adaptor to in the registry.
  Arg [3]    : The adaptor to be added to the registry.
  Example    : Bio::EnsEMBL::Registry->add_db($db, "lite", $dba);
  Returntype : none
  Exceptions : none
  Status     : At Risk.
             : This is here for backwards compatibility only and may
             : be removed eventually.  Solution is to make sure the
             : db and the adaptor have the same species and the call
             : is then no longer needed.

=cut

sub add_db {
  my ( $self, $db, $name, $adap ) = @_;
  my $storage = $self->storage();
  #No warnings brought in due to some overzelous webcode triggering a lot of warnings.
  no warnings 'uninitialized';
  my $lc_db_species = lc($db->species());
  if ( $lc_db_species ne lc( $adap->species ) ) {
    my $lc_group = lc($db->group());
    my $lc_name = lc($name);
    $storage->{_SPECIES}{ $lc_db_species }{ $lc_group }{'_special'}{ $lc_name } = $adap;
  }
  return;
}

=head2 remove_db

  Arg [1]    : db (DBAdaptor) to remove adaptor from.
  Arg [2]    : name to remove the adaptor from in the registry.
  Example    : my $db = Bio::EnsEMBL::Registry->remove_db($db, "lite");
  Returntype : adaptor The instance which has just been removed from the registry
  Exceptions : none
  Status     : At Risk.
             : This is here for backwards compatibility only and may
             : be removed eventually.  Solution is to make sure the
             : db and the adaptor have the same species and the call
             : is then no longer needed.

=cut

sub remove_db {
  my ( $self, $db, $name ) = @_;
  my $storage = $self->storage();
  my $lc_db_species = lc($db->species());
  my $lc_group = lc($db->group());
  my $lc_name = lc($name);
  my $ret = $storage->{_SPECIES}{ $lc_db_species }{ $lc_group }{'_special'}{ $lc_name };
  $storage->{_SPECIES}{ $lc_db_species }{ $lc_group }{'_special'}{ $lc_name } = undef;
  return $ret;
}

=head2 get_db

  Arg [1]    : db (DBAdaptor) to get adaptor from.
  Arg [2]    : name to get the adaptor for in the registry.
  Example    : my $db = Bio::EnsEMBL::Registry->get_db("Human", "core", "lite");
  Returntype : adaptor
  Exceptions : See get_DBAdaptor()
  Status     : At Risk.
             : This is here for backwards compatibility only and may
             : be removed eventually.  Solution is to make sure the
             : db and the adaptor have the same species then call
             : get_DBAdaptor instead.

=cut

sub get_db {
  my ( $self, $db, $name ) = @_;
  my $storage = $self->storage();
  my $lc_species = lc($db->species());
  my $lc_group = lc($db->group());
  my $ret = Bio::EnsEMBL::Registry->get_DBAdaptor($lc_species, $lc_group);
  return $ret if defined $ret;
  my $lc_name = lc($name);
  return $storage->{_SPECIES}{ $lc_species }{ $lc_group }{'_special'}{ $lc_name };
}

=head2 get_all_db_adaptors

  Arg [1]    : db (DBAdaptor) to get all the adaptors from.
  Example    : my $db = Bio::EnsEMBL::Registry->get_all_db_adaptors($db);
  Returntype : adaptor
  Exceptions : none
  Status     : At Risk.
             : This is here for backwards compatibility only and
             : may be removed eventually.  Solution is to make
             : sure the dbs all have the same species then call
             : get_all_DBAdaptors(-species => "human");


=cut

sub get_all_db_adaptors {
  my ( $self, $db ) = @_;
  my %ret = ();
  my $storage = $self->storage();
  my $lc_db_species = lc($db->species());
  my $lc_group = lc($db->group());
  
  # we now also want to add all the DBAdaptors for the same species.
  # as add_db_adaptor does not add if it is from the same species.

  foreach my $dba ( @{ $storage->{'_DBA'} } ) {
    if ( lc( $dba->species() ) eq $lc_db_species ) {
      $ret{ $dba->group() } = $dba;
    }
  }

  my $alias = $self->get_alias($lc_db_species);
  foreach my $key (keys %{ $storage->{_SPECIES}{ $alias }{ $lc_group }{'_special'} }) {
    $ret{$key} = $storage->{_SPECIES}{ $alias }{ $lc_group }{'_special'}{ $key };
  }

  return \%ret;
} ## end sub get_all_db_adaptors


#
# DBAdaptors
#

=head2 add_DBAdaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : DBAdaptor to be added to the registry.
  Example    : Bio::EnsEMBL::Registry->add_DBAdaptor("Human", "core", $dba);
  Returntype : none
  Exceptions : none
  caller     : internal
  Status     : Stable

=cut

sub add_DBAdaptor {
  my ( $self, $species, $group, $adap ) = @_;
  my $storage = $self->storage();
  if ( !( $self->alias_exists($species) ) ) {
    $self->add_alias( $species, $species );
  }
  $species = $self->get_alias($species);
  $storage->{_SPECIES}{$species}{ lc($group) }{'_DB'} = $adap;
  if ( !defined( $storage->{'_DBA'} ) ) {
    $storage->{'_DBA'} = [$adap];
  } 
  else {
    push( @{ $storage->{'_DBA'} }, $adap );
  }
  return;
}

=head2 get_DBAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Arg [3]    : if set will not give warnings when looking for alias.
  Example    : $dba = Bio::EnsEMBL::Registry->get_DBAdaptor("Human", "core");
  Returntype : DBAdaptor
  Exceptions : If $species is not defined and if no valid internal name 
               could be found for $species. If thrown check your API and DB
               version 
  Status     : Stable

=cut

sub get_DBAdaptor {
  my ( $self, $species, $group, $no_alias_check ) = @_;

  if ( !defined($species) ) {
    throw('Species not defined.');
  }

  my $ispecies = $self->get_alias( $species, $no_alias_check );

  if ( !defined($ispecies) ) {
    if(! $no_alias_check) {
      throw("Can not find internal name for species '$species'");
    }
  }
  else { $species = $ispecies }
  my $storage = $self->storage();
  return $storage->{_SPECIES}{$species}{ lc($group) }{'_DB'};
}

=head2 get_all_DBAdaptors

  Arg [SPECIES]: (optional) string 
                  species name to get adaptors for
  Arg [GROUP]  : (optional) string 
                  group name to get adaptors for
  Example      : 
                @dba =
                  @{ Bio::EnsEMBL::Registry->get_all_DBAdaptors() };

                @human_dbas =
                  @{ Bio::EnsEMBL::Registry->get_all_DBAdaptors(
                    -species => 'human'
                  ) };

  Returntype   : list of DBAdaptors
  Exceptions   : none
  Status       : Stable

=cut

sub get_all_DBAdaptors {
  my ( $self, @args ) = @_;

  my ( $species, $group ) = rearrange( [qw(SPECIES GROUP)], @args );

  if ( defined($species) ) { $species = $self->get_alias($species) }
  if ( defined($group) ) { $group = lc($group) }

  my @ret;
  my $storage = $self->storage();
  foreach my $dba ( @{ $storage->{'_DBA'} } ) {
    if ( ( !defined($species) || $species eq lc( $dba->species() ) )
      && ( !defined($group) || $group eq lc( $dba->group() ) ) )
    {
      push( @ret, $dba );
    }
  }

  return \@ret;
}

=head2 get_all_DBAdaptors_by_connection

  Arg [1]    : DBConnection used to find DBAdaptors
  Returntype : reference to list of DBAdaptors
  Exceptions : none
  Example    : @dba = @{ Bio::EnsEMBL::Registry
                  ->get_all_DBAdaptors_by_connection($dbc) };
  Status     : Stable

=cut

sub get_all_DBAdaptors_by_connection {
  my ( $self, $dbc_orig ) = @_;
  my @return;
  my $storage = $self->storage();
  foreach my $dba ( @{ $storage->{'_DBA'} } ) {
    my $dbc = $dba->dbc();
    if (defined $dbc && $dbc->equals($dbc_orig) ) {
      push( @return, $dba );
    }
  }
  return \@return;
}

=head2 get_all_DBAdaptors_by_dbname

  Arg [1]    : string, name of database
  Returntype : reference to list of DBAdaptors
  Exceptions : none
  Example    : @dba = @{ Bio::EnsEMBL::Registry
                  ->get_all_DBAdaptors_by_dbname($dbname) };
  Status     : Stable

=cut

sub get_all_DBAdaptors_by_dbname {
  my ( $self, $dbname ) = @_;
  my @return;
  my $storage = $self->storage();
  foreach my $dba ( @{ $storage->{'_DBA'} } ) {
    my $dbc = $dba->dbc();
    if ( defined($dbc) && $dbc->dbname() eq $dbname ) {
      push( @return, $dba );
    }
  }
  return \@return;
}

=head2 remove_DBAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Example    : $dba = Bio::EnsEMBL::Registry->remove_DBAdaptor("Human", "core");
  Returntype : none
  Exceptions : none
  Status     : At risk

=cut

sub remove_DBAdaptor {
  my ( $self, $species, $group ) = @_;
  $species = $self->get_alias($species);
  my $storage = $self->storage();
  delete $storage->{_SPECIES}{$species}{$group};
  # This will remove the DBAdaptor and all the other adaptors
  # Now remove if from the _DBA array
  my $index;
  foreach my $i ( 0 .. $#{ $storage->{'_DBA'} } ) {
    my $dba = $storage->{'_DBA'}->[$i];
    if ( $dba->species() eq $species && $dba->group() eq $group ) {
      $index = $i;
      last;
    }
  }

  # Now remove from _DBA cache
  if ( defined($index) ) {
    splice( @{ $storage->{'_DBA'} }, $index, 1 );
  }

  return;
} ## end sub remove_DBAdaptor


=head2 reset_DBAdaptor

  Arg [1]     : String - species e.g. homo_sapiens
  Arg [2]     : String - DB group e.g. core
  Arg [3]     : String - new dbname
  Arg [4]     : (optional) String - new host to use
  Arg [5]     : (optional) Integer - new port to use
  Arg [6]     : (optional) String - new username to use
  Arg [7]     : (optional) String - new password to use
  Arg [8]     : HashRef - Provide additional key-value parameters to pass into the 
                DBAdaptor's new constructor e.g. eFG dnadb params for auto selecting dnadb
  Usage       : $reg->reset_registry_db( 'homo_sapiens', 'core', 'homo_sapiens_core_37_35j' );
  Description : Resets a DB within the registry.
  Exceptions  : Throws if mandatory params not supplied
                Throws if species name is not already seen by the registry
                Throws if no current DB for species/group available
  Status      : At risk

=cut

sub reset_DBAdaptor {
  my (
    $self, $species, $group, $dbname, $host,
    $port, $user,    $pass,  $params
  ) = @_;

  # Check mandatory params
  throw "Must provide a species to redefine a DB in the registry" if ! defined $species;
  throw "Must provide a group to redefine a DB in the registry" if ! defined $group;
  throw "Must provide a dbname to redefine a DB in the registry" if ! defined $dbname;

  # Validate species here
  my $alias = $self->get_alias($species);
  throw("Could not find registry alias for species:\t$species")
    if ( !defined $alias );

  # Get all current defaults if not defined
  my $db = $self->get_DBAdaptor( $alias, $group );
  my $class;

  if ($db) {
    $class = ref($db);
    my $dbc = $db->dbc();
    $host ||= $dbc->host;
    $port ||= $dbc->port;
    $user ||= $dbc->username;
    $pass ||= $dbc->password;
  } 
  #Is this really re-defining? If the DBAdaptor isn't there then this 
  #is just construction by another name
  else {
    #Now we need to test mandatory params
    throw "No comparable $alias $group DB present in the registry. You must provide a dbhost" unless $host;
    throw "No comparable $alias $group DB present in the registry. You must provide a dbuser" unless $user;
    $class = $group2adaptor{ lc($group) };
  }
  $self->remove_DBAdaptor( $alias, $group );

  # ConfigRegistry should automatically add this to the Registry
  $db = $class->new(
    -user    => $user,
    -host    => $host,
    -port    => $port,
    -pass    => $pass,
    -dbname  => $dbname,
    -species => $alias,
    -group   => $group,
    %{$params} );

  return $db;
} ## end sub reset_DBAdaptor


#
# DNA Adaptors
#

=head2 add_DNAAdaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the species to get the dna from
  Arg [4]    : name of the group to get the dna from
  Example    : Bio::EnsEMBL::Registry->add_DNAAdaptor("Human", "estgene", "Human", "core");
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub add_DNAAdaptor {
  my ( $self, $species, $group, $dnadb_species, $dnadb_group ) = @_;

  $species       = $self->get_alias($species);
  $dnadb_species = $self->get_alias($dnadb_species);
  if ( $dnadb_group->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') ) {
    deprecated("");
  } else {
    my $storage = $self->storage();
    $storage->{_SPECIES}{$species}{ lc($group) }{'_DNA'} = $dnadb_group;
    $storage->{_SPECIES}{$species}{ lc($group) }{'_DNA2'} = $dnadb_species;
  }
  return;
}

=head2 get_DNAAdaptor

  Arg [1]    : name of the species to get the adaptor for in the registry.
  Arg [2]    : name of the group to get the adaptor for in the registry.
  Example    : $dnaAdap = Bio::EnsEMBL::Registry->get_DNAAdaptor("Human", "core");
  Returntype : adaptor
  Exceptions : none
  Status     : Stable

=cut

sub get_DNAAdaptor {
  my ( $self, $species, $group ) = @_;

  $species = $self->get_alias($species);
  my $lc_group = lc($group);
  my $storage = $self->storage();
  my $new_group = $storage->{_SPECIES}{$species}{ $lc_group }{'_DNA'};
  my $new_species = $storage->{_SPECIES}{$species}{ $lc_group }{'_DNA2'};

  if ( defined $new_group ) {
    return $self->get_DBAdaptor( $new_species, $new_group );
  }

  return;
}

#
# General Adaptors
#

=head2 add_adaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the type to add the adaptor to in the registry.
  Arg [4]    : The DBAdaptor to be added to the registry.
  Arg [5]    : (optional) Set to allow overwrites of existing adaptors.
  Example    : Bio::EnsEMBL::Registry->add_adaptor("Human", "core", "Gene", $adap);
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub add_adaptor {
  my ( $self, $species, $group, $type, $adap, $reset ) = @_;

  $species = $self->get_alias($species);
  my $lc_group = lc($group);
  my $lc_type = lc($type);

  my $storage = $self->storage();

  # Since the adaptors are not stored initially, only their class paths
  # when the adaptors are obtained, we need to store these instead.  It
  # is not necessarily an error if the registry is overwritten without
  # the reset set but it is an indication that we are overwriting a
  # database which should be a warning for now
  
  if ( defined($reset) )
  {    # JUST RESET THE HASH VALUE NO MORE PROCESSING NEEDED
    $storage->{_SPECIES}{$species}{ $lc_group }{ $lc_type } =
      $adap;
    return;
  }

  if (
    defined(
      $storage->{_SPECIES}{$species}{ $lc_group }{ $lc_type }
    ) )
  {
  # print STDERR (
  #      "Overwriting Adaptor in Registry for $species $group $type\n");
    $storage->{_SPECIES}{$species}{ $lc_group }{ $lc_type } =
      $adap;
    return;
  }
  $storage->{_SPECIES}{$species}{ $lc_group }{ $lc_type } =
    $adap;

  if ( !defined( $storage->{_SPECIES}{$species}{'list'} ) ) {
    $storage->{_SPECIES}{$species}{'list'} = [$type];
  } else {
    push( @{ $storage->{_SPECIES}{$species}{'list'} }, $type );
  }

  if ( !defined( $storage->{_TYPE}{ $lc_type }{$species} ) ) {
    $storage->{_TYPE}{ $lc_type }{$species} = [$adap];
  } else {
    push( @{ $storage->{_TYPE}{ $lc_type }{$species} },
      $adap );
  }
  return;
} ## end sub add_adaptor

sub add_all_adaptors {
  my ($self, $dba) = @_;

  my $species = $self->get_alias($dba->species());
  my ($group) = lc($dba->group());
  
  my $storage = $self->storage();
  my $adaptors_hash = $dba->get_available_adaptors(); 
  my @types = keys %{$adaptors_hash};
    
  foreach my $type (@types) {
    my $lc_type = lc($type);
    my $adaptor = $adaptors_hash->{$type};
    
    #Add the adaptor into the registry
    $storage->{_SPECIES}->{$species}->{$group}->{$lc_type} = $adaptor;
    
    #Add adaptor to the types hash
    if(!exists $storage->{_TYPE}->{$lc_type}->{$species}) {
      $storage->{_TYPE}->{$lc_type}->{$species} = [];
    }
    push(@{$storage->{_TYPE}->{$lc_type}->{$species}}, $adaptor);
  }
  
  # Add a the list of types the _SPECIES hash
  $storage->{_SPECIES}{$species}{'list'} = \@types;
  
  return;
}

=head2 get_adaptor

  Arg [1]    : name of the species to add the adaptor to in the registry.
  Arg [2]    : name of the group to add the adaptor to in the registry.
  Arg [3]    : name of the type to add the adaptor to in the registry.
  Example    : $adap = Bio::EnsEMBL::Registry->get_adaptor("Human", "core", "Gene");
  Returntype : adaptor
  Exceptions : Thrown if a valid internal name cannot be found for the given 
               name. If thrown check your API and DB version. Also thrown if
               no type or group was given
  Status     : Stable

=cut

sub get_adaptor {
  my ( $self, $species, $group, $type ) = @_;
  
  my $ispecies = $self->get_alias($species);

  if ( !defined($ispecies) ) {
    throw("Can not find internal name for species '$species'");
  }
  else { $species = $ispecies }
  
  throw 'No adaptor group given' if ! defined $group;
  throw 'No adaptor type given' if ! defined $type;
  
  my $lc_group = lc($group);
  my $lc_type = lc($type);
  
  if($lc_type =~ /adaptor$/) {
    warning("Detected additional Adaptor string in given the type '$type'. Removing it to avoid possible issues. Alter your type to stop this message");
    $type =~ s/adaptor$//i;
    $lc_type = lc($type);
  }

  # For historical reasons, allow use of group 'regulation' to refer to
  # group 'funcgen'.
  if ( $lc_group eq 'regulation' ) { $lc_group = 'funcgen' }

  my %dnadb_adaptors = (
    'sequence'                 => 1,
    'assemblymapper'           => 1,
    'karyotypeband'            => 1,
    'repeatfeature'            => 1,
    'coordsystem'              => (($lc_group ne 'funcgen') ? 1 : undef),
    'assemblyexceptionfeature' => 1
  );
  
  my $storage = $self->storage();
  
  my $dnadb_group =
    $storage->{_SPECIES}{$species}{ $lc_group }{'_DNA'};

  if ( defined($dnadb_group)
    && defined( $dnadb_adaptors{ $lc_type } ) )
  {
    $species =
      $storage->{_SPECIES}{$species}{ $lc_group }{'_DNA2'};
    $group = $dnadb_group;
  }

  my $ret =
    $storage->{_SPECIES}{$species}{ $lc_group }{ $lc_type };

  if ( !defined($ret) ) { return }
  if ( ref($ret) )      { return $ret }

  # Not instantiated yet

  my $dba = $storage->{_SPECIES}{$species}{ $lc_group }{'_DB'};
  my $module = $ret;

  my $test_eval = eval "require $module"; ## no critic
  if ($@ or (!$test_eval)) {
    warning("'$module' cannot be found.\nException $@\n");
    return;
  }

  if (
    !defined(
      $storage->{_SPECIES}{$species}{ $lc_group }{'CHECKED'} )
    )
  {
    $storage->{_SPECIES}{$species}{ $lc_group }{'CHECKED'} = 1;
    $self->version_check($dba);
  }

  my $adap = "$module"->new($dba);
  Bio::EnsEMBL::Registry->add_adaptor( $species, $lc_group, $type, $adap,
                                       'reset' );
  $ret = $adap;

  return $ret;
} ## end sub get_adaptor

=head2 get_all_adaptors

  Arg [SPECIES] : (optional) string 
                  species name to get adaptors for
  Arg [GROUP] : (optional) string 
                  group name to get adaptors for
  Arg [TYPE] : (optional) string 
                  type to get adaptors for
  Example    : @adaps = @{Bio::EnsEMBL::Registry->get_all_adaptors()};
  Returntype : ref to list of adaptors
  Exceptions : none
  Status     : Stable

=cut

sub get_all_adaptors{
  my ($self,@args)= @_;
  my ($species, $group, $type);
  my @ret=();
  my (%species_hash, %group_hash, %type_hash);


  if(@args == 1){ # Old species only one parameter
    warn("-SPECIES argument should now be used to get species adaptors");
    $species = $args[0];
  }
  else{
    # new style -SPECIES, -GROUP, -TYPE
    ($species, $group, $type) =
      rearrange([qw(SPECIES GROUP TYPE)], @args);
  }
  
  my $storage = $self->storage();

  if(defined($species)){
    $species_hash{$species} = 1;
  }
  else{
    # get list of species
    foreach my $dba (@{$storage->{'_DBA'}}){
      $species_hash{lc($dba->species())} = 1;
    }
  }
  if(defined($group)){
    $group_hash{$group} = 1;
  }
  else{
    foreach my $dba (@{$storage->{'_DBA'}}){
      $group_hash{lc($dba->group())} = 1;
    }
  }

  if ( defined($type) ) {
    $type_hash{$type} = 1;
  } else {
    foreach my $dba ( @{ $storage->{'_DBA'} } ) {
      foreach my $ty (
        @{ $storage->{_SPECIES}{ lc( $dba->species ) }{'list'} }
        )
      {
        $type_hash{ lc($ty) } = 1;
      }
    }
  }

  ### NOW NEED TO INSTANTIATE BY CALLING get_adaptor
  foreach my $sp ( keys %species_hash ) {
    foreach my $gr ( keys %group_hash ) {
      foreach my $ty ( keys %type_hash ) {
        my $temp = $self->get_adaptor( $sp, $gr, $ty );
        if ( defined($temp) ) {
          push @ret, $temp;
        }
      }
    }
  }

  return (\@ret);
}

=pod2 find_free_species_name

  Arg[1]      : String
  Description : Scans through the registry looking for 
  Example     : my $name = $registry->find_free_species_name('species'); # returns species
                $registry->add_alias($name);
                $name = $registry->find_free_species_name('species'); # returns species1 because species has been added
=cut

sub find_free_species_name {
  my ($self, $species) = @_;
  my $storage = $self->storage();
  my $i = 0;
  while(1) {
    if($storage->alias_exists($species)) {
      $i++;
      next;
    }
    last;
  }
  return ($i == 0) ? $species : $species.$i;
}

=head2 add_alias

  Arg [1]    : name of the species to add alias for
  Arg [2]    : String or ArrayRef name (or names) of the alias(es) to add
  Example    : Bio::EnsEMBL::Registry->add_alias("Homo Sapiens","Human");
               Bio::EnsEMBL::Registry->add_alias("Homo Sapiens",['human', 'homo_sapiens', '9606']);
  Description: Add alternative names for a species.
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub add_alias {
  my ($self, $species, $key) = @_;
  my $lc_species = lc($species);
  my $storage = $self->storage();
  if(defined $key && ref($key) eq 'ARRAY') {
    foreach my $alias (map {lc($_)} @{$key}) {
      $storage->{'_ALIAS'}{$alias} = $lc_species;
    }
  }
  else {
    $storage->{'_ALIAS'}{lc($key)} = $lc_species;
  }
  return;
}

=head2 remove_alias

  Arg [1]    : name of the species to remove alias for
  Arg [2]    : name of the alias
  Example    : Bio::EnsEMBL::Registry->remove_alias("Homo Sapiens","Human");
  Description: remove alternative name for the species.
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub remove_alias{
  my ($self, $species, $key) = @_;
  delete $self->storage->{'_ALIAS'}{lc($key)};
  return;
}

=head2 get_alias

  Arg [1]    : name of the possible alias to get species for
  Example    : Bio::EnsEMBL::Registry->get_alias("Human");
  Description: get proper species name.
  Returntype : species name
  Exceptions : none
  Status     : Stable

=cut

sub get_alias {
  my ( $self, $key, $no_warn ) = @_;
  my $lc_key = lc($key);
  my $storage = $self->storage();
  if ( !defined( $storage->{'_ALIAS'}{ $lc_key } ) ) {
    if ( ( !defined( $storage->{_SPECIES}{ $lc_key } ) ) and
         ( !defined( $storage->{_ALIAS}{ $lc_key } ) ) )
    {
      if ( ( !defined($no_warn) ) or ( !$no_warn ) ) {
        warning( "$key is not a valid species name " .
                 "(check DB and API version)" );
      }
      return;
    }
    else { return $key }
  }

  return $storage->{'_ALIAS'}{ $lc_key };
}

=head2 get_all_aliases

  Arg [1]    : Species name to retrieve aliases for
               (may be an alias as well).
  Example    : Bio::EnsEMBL::Registry->get_all_aliases('Homo sapiens');
  Description: Returns all known aliases for a given species (but not the
               species name/alias that was given).
  Returntype : ArrayRef of all known aliases
  Exceptions : none
  Status     : Development

=cut

sub get_all_aliases {
  my ( $self, $key ) = @_;
  my $lc_key = lc($key);
  my $storage = $self->storage();
  my $species = $storage->{_ALIAS}{ $lc_key };

  my @aliases;
  if ( defined($species) ) {
    foreach my $alias ( keys( %{ $storage->{_ALIAS} } ) ) {
      if ( $species ne $alias
        && $species eq $storage->{_ALIAS}{ $alias } )
      {
        push( @aliases, $alias );
      }
    }
  }

  return \@aliases;
}

=head2 alias_exists

  Arg [1]    : name of the possible alias to get species for
  Example    : Bio::EnsEMBL::Registry->alias_exists("Human");
  Description: does the species name exist.
  Returntype : 1 if exists else 0
  Exceptions : none
  Status     : Stable

=cut

sub alias_exists {
  my ( $self, $key ) = @_;
  return ( defined $self->storage()->{'_ALIAS'}{ lc($key) } ) ? 1 : 0;
}

=head2 set_disconnect_when_inactive

  Example    : Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
  Description: Set the flag to make sure that the database connection is dropped if
               not being used on each database.
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub set_disconnect_when_inactive{
  my ($self) = @_;
  foreach my $dba ( @{$self->get_all_DBAdaptors()}){
    my $dbc = $dba->dbc;
    # Disconnect if connected
    $dbc->disconnect_if_idle() if $dbc->connected();
    $dbc->disconnect_when_inactive(1);
  }
  return;
}

=head2 set_reconnect_when_lost

  Example    : Bio::EnsEMBL::Registry->set_reconnect_when_lost();
  Description: Set the flag to make sure that the database connection is not lost before it is used.
               This is useful for long running jobs (over 8hrs).
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub set_reconnect_when_lost{
  my ($self) = @_;
  foreach my $dba ( @{$self->get_all_DBAdaptors()}){
    my $dbc = $dba->dbc;
    $dbc->reconnect_when_lost(1);
  }
  return;
}

=head2 disconnect_all

  Example    : Bio::EnsEMBL::Registry->disconnect_all();
  Description: disconnect from all the databases.
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub disconnect_all {
  my ($self) = @_;
  foreach my $dba ( @{$self->get_all_DBAdaptors()} ){
    my $dbc = $dba->dbc;
    next unless $dbc;
    # Disconnect if connected
    $dbc->disconnect_if_idle() if $dbc->connected();
  }
  return;
}

=head get_DBAdaptor_count

  Example     : Bio::EnsEMBL::Registry->get_DBAdaptor_count();
  Description : Returns the count of database adaptors currently held by 
                the registry
  Returntype  : Int count of database adaptors currently known
  Exceptions  : None
 
=cut

sub get_DBAdaptor_count {
  my ($self) = @_;
  my $storage = $self->storage();
  return scalar(@{$storage->{'_DBA'}}) if(defined $storage->{'_DBA'});
  return 0;
}

=head2 change_access

  Will change the username and password for a set of databases.
  if host,user or database names are missing then these are not checked.
  So for example if you do not specify a database then ALL databases on
  the specified  host and port will be changed.

  Arg [1]    : name of the host to change access on
  Arg [2]    : port number to change access on
  Arg [3]    : name of the user to change access on
  Arg [4]    : name of the database to change access on
  Arg [5]    : name of the new user
  Arg [6]    : new password

  Example    : Bio::EnsEMBL::Registry->get_alias("Human");
  Description: change username and password on one or more databases
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut

sub change_access{
  my ($self, $host,$port,$user,$dbname,$new_user,$new_pass) = @_;
  my $storage = $self->storage();
  foreach my $dba ( @{$storage->{'_DBA'}}){
    my $dbc = $dba->dbc;
    if((((!defined($host)) or ($host eq $dbc->host))) and
       (((!defined($port)) or ($port eq $dbc->port))) and
       (((!defined($user)) or ($user eq $dbc->username))) and
       ((!defined($dbname)) or ($dbname eq $dbc->dbname))){
      if($dbc->connected()){
        $dbc->db_handle->disconnect();
        $dbc->connected(undef);
      }
      # over write the username and password
      $dbc->username($new_user);
      $dbc->password($new_pass);
    }
  }
  return;
}

=head2 load_registry_from_url

  Arg [1] : string $url
  Arg [2] : (optional) integer
            If not 0, will print out all information.
  Arg [3] : (optional) integer
          This option will turn off caching for slice features, so,
          every time a set of features is retrieved, they will come
          from the database instead of the cache. This option is only
          recommended for advanced users, specially if you need to
          store and retrieve features. It might reduce performance when
          querying the database if not used properly. If in doubt, do
          not use it or ask in the developer mailing list.

  Example : load_registry_from_url(
            'mysql://anonymous@ensembldb.ensembl.org:3306');
            
            load_registry_from_url(
            'mysql://anonymous@ensembldb.ensembl.org:3306/homo_sapiens_core_65_37?group=core&species=homo_sapiens'
            );
            
            load_registry_from_url(
            'mysql://anonymous@ensembldb.ensembl.org:3306/homo_sapiens_core_65_37?group=core'
            );
            

  Description: Will load the correct versions of the ensembl
               databases for the software release it can find on
               a database instance into the registry. Also adds
               a set of standard aliases. The url format is:
               mysql://[[username][:password]@]hostname[:port].  You
               can also request a specific version for the databases
               by adding a slash and the version number but your
               script may crash as the API version won't match the
               DB version.
               
               You can also specify a database name which will cause the 
               loading of a single DBAdaptor instance. Parameters are
               mapped from a normal URL parameter set to their DBAdaptor
               equivalent. Group must be defined.
               
  Returntype : Int count of the DBAdaptor instances which can be found in the 
               registry

  Exceptions : Thrown if the given URL does not parse according to the above 
               scheme and if the specified database cannot be connected to 
               (see L<load_registry_from_db> for more information)
  Status     : Stable
 
=cut

sub load_registry_from_url {
  my ( $self, $url, $verbose, $no_cache ) = @_;
  my $original_count = $self->get_DBAdaptor_count();
  my $loader = Bio::EnsEMBL::Registry::Loader->new($self, $verbose, $no_cache);
  $loader->load_registry($url);
  my $count = $self->get_DBAdaptor_count() - $original_count;
  return $count >= 0 ? $count : 0; 
}

=head2 load_registry_from_db

  Arg [HOST] : string
                The domain name of the database host to connect to.

  Arg [USER] : string
                The name of the database user to connect with.

  Arg [PASS] : (optional) string
                The password to be used to connect to the database.

  Arg [PORT] : (optional) integer
                The port to use when connecting to the database.

  Arg [VERBOSE]: (optional) boolean
                Whether to print database messages. This includes a listing
                of all available species & databases.

  Arg [SPECIES]: (optional) string
                By default, all databases that are found on the
                server and that corresponds to the correct release
                are probed for aliases etc.  For some people,
                depending on where they are in the world, this might
                be a slow operation.  With the '-species' argument,
                one may reduce the startup time by restricting the
                set of databases that are probed to those of a
                particular species.

                Note that the latin name of the species is required,
                e.g., 'homo sapiens', 'gallus gallus', 'callithrix
                jacchus' etc.  It may be the whole species name,
                or only the first part of the name, e.g. 'homo',
                'gallus', or 'callithrix'.  This will be used in
                matching against the name of the databases.

  Arg [DB_VERSION]: (optional) integer
                By default, only databases corresponding to the
                current API version are loaded.  This argument
                allows the script to use databases from another
                version although it might not work properly.  This
                argument should only be used for production or
                testing purposes and if you really know what you are
                doing.

  Arg [WAIT_TIMEOUT]: (optional) integer
                Time in seconds for the wait timeout to happen.
                Time after which the connection is deleted if not
                used.  By default this is 28800 (8 hours), so set
                this to greater than this if your connection are
                getting deleted.  Only set this if you are having
                problems and know what you are doing.

   Arg [-NO_CACHE]: (optional) boolean
                This option will turn off caching for slice features,
                so, every time a set of features is retrieved, they
                will come from the database instead of the cache.  This
                option is only recommended for advanced users, specially
                if you need to store and retrieve features.  It might
                reduce performance when querying the database if not
                used properly.  If in doubt, do not use it or ask in the
                developer mailing list.

   Arg [SPECIES_SUFFIX]: (optional) string
                This option will append the string to the species name
                in the registry for all databases found on this server.

  Example :

    $registry->load_registry_from_db(
      -host    => 'ensembldb.ensembl.org',
      -user    => 'anonymous',
      -verbose => '1'
    );

  Description: Will load the correct versions of the Ensembl
               databases for the software release it can find on a
               database instance into the registry.  Also adds a set
               of standard aliases.

  Returntype : Int count of the DBAdaptor instances which can be found in the 
               registry due to this method call.

  Exceptions : Thrown if the given MySQL database cannot be connected to
               or there is any error whilst querying the database.
  Status     : Stable

=cut

sub load_registry_from_db {
  my ( $self, @args ) = @_;
  my $original_count = $self->get_DBAdaptor_count();
  my $loader = Bio::EnsEMBL::Registry::Loader::Database->new($self);
  $loader->load_registry(@args);
  my $count = $self->get_DBAdaptor_count() - $original_count;
  return $count >= 0 ? $count : 0; 
}

=head2 find_and_add_aliases

  Arg [ADAPTOR] : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor
                  The adaptor to use to retrieve aliases from.

  Arg [GROUP]   : (optional) string
                  The group you want to find aliases for. If not
                  given assumes all types.

  Arg [HANDLE]  : (optional) DBI database handle
                  A connected database handle to use instead of
                  the database handles stored in the DBAdaptors.
                  Bypasses the use of MetaContainer.

  Arg [SPECIES_SUFFIX]: (optional) string
                  This option will append the string to the species
                  name in the registry for all databases.

  Example       : Bio::EnsEMBL::Registry->find_and_add_aliases(
                    -ADAPTOR => $dba,
                    -GROUP   => 'core'
                  );

  Description   : Looks in the meta container for each database for
                  an entry called "species.alias".  If any are found
                  then the species adaptor is registered to that
                  set of aliases.  This can work across any adaptor
                  which has a MetaContainer.  If no MetaContainer
                  can be returned from a given adaptor then no alias
                  searching is performed.

  Return type   : none
  Exceptions    : Throws if an alias is found in more than one species.
  Status        : Stable

=cut

sub find_and_add_aliases {
  my $self = shift @_;
  my ($adaptor, $group, $dbh, $species_suffix ) = 
    rearrange( [ 'ADAPTOR', 'GROUP', 'HANDLE', 'SPECIES_SUFFIX' ], @_ );

  $species_suffix ||=  q{};

  my @dbas;
  if ( defined($adaptor) ) {
    @dbas = ($adaptor);
  } 
  elsif ( defined($dbh) ) {
    if($species_suffix) {
      @dbas = grep { $_->species() =~ /$species_suffix/ } @{$self->get_all_DBAdaptors(-GROUP => $group)};
    }
    else {
      @dbas = @{ $self->get_all_DBAdaptors( '-GROUP' => $group ) };
    }
  }
  else {
    @dbas = @{ $self->get_all_DBAdaptors( '-GROUP' => $group ) };
  }
  
  my $aliases_for_dbc = {};

  foreach my $dba (@dbas) {
    my @aliases;
    my $species = $dba->species();

    if ( defined($dbh) ) {
    	
    	my $dbname = $dba->dbc()->dbname();

          if (!defined $aliases_for_dbc->{$dbname}) {

                my $sth = $dbh->prepare(sprintf("SELECT species_id,meta_value FROM %s.meta " 
                    . "WHERE meta_key = 'species.alias' ", $dbh->quote_identifier($dbname))
                );

                # Execute, and don't care about errors (there will be errors for
                # databases without a 'meta' table.
                $sth->{'PrintError'} = 0;
                $sth->{'RaiseError'} = 0;
                if (!$sth->execute()) { next }
                $sth->{'PrintError'} = $dbh->{'PrintError'};
                $sth->{'RaiseError'} = $dbh->{'RaiseError'};

                my $alias;
                my $species_id;
                $sth->bind_columns(\$species_id, \$alias);
                while ($sth->fetch()) {
                  push(@{$aliases_for_dbc->{$dbname}{$species_id}}, $alias);
                }
          }

          @aliases = @{$aliases_for_dbc->{$dbname}{$dba->species_id()}||[]}

    } else {
      my $meta_container = eval { $dba->get_MetaContainer() };

      if ( defined($meta_container) ) {
        push( @aliases,
              @{ $meta_container->list_value_by_key('species.alias') }
        );
      }

      # Need to disconnect so we do not spam the MySQL servers trying to
      # get aliases.  Can only call disonnect if dbc was defined.
      if ( defined( $dba->dbc() ) ) {
        $dba->dbc()->disconnect_if_idle();
      }
    }

    foreach my $alias (@aliases) {
      my $alias_suffix = $alias.$species_suffix;
      #Lowercase because stored aliases are lowercased
      my $lc_species = lc($species);
      my $lc_alias_suffix = lc($alias_suffix);
      if (   !$self->alias_exists( $alias_suffix )
           && $lc_species ne $lc_alias_suffix )
      {
        $self->add_alias( $species, $alias_suffix );
      } elsif (
             $lc_species ne $self->get_alias( $alias_suffix ) )
      {
        $self->remove_alias( $species, $alias_suffix );
      }
    }

  } ## end foreach my $dba (@dbas)
  return;
} ## end sub find_and_add_aliases


=head2 load_registry_from_multiple_dbs

  Arg [1]   : Array of hashes, each hash being a set of arguments to
              load_registry_from_db() (see above).

  Example   :

    $registry->load_registry_from_multiple_dbs( {
        '-host'    => 'ensembldb.ensembl.org',
        '-user'    => 'anonymous',
        '-verbose' => '1'
      },
      {
        '-host'     => 'server.example.com',
        '-user'     => 'anonymouse',
        '-password' => 'cheese',
        '-verbose'  => '1'
      } );

  Description:  Will call load_registry_from_db() (see above)
                multiple times and merge the resulting registries
                into one, effectively allowing a user to connect to
                databases on multiple database servers from within
                one program.

                If a database is found on more than one server, the
                first found instance of that database will be used.

  Returntype : Int count of the DBAdaptor instances which can be found in the 
               registry

=cut

sub load_registry_from_multiple_dbs {
  my ( $self, @args ) = @_;
  
  my $current_storage = $self->storage();

  my $original_count = $self->get_DBAdaptor_count();

  my %merged_register = %registry_register;

  foreach my $arg (@args) {
    local %registry_register = ();

    my $verbose;

    ($verbose) = rearrange( ['VERBOSE'], %{$arg} );

    $self->load_registry_from_db( %{$arg} );

    #
    # Merge the localized %registry_register into %merged_register.
    #

    # Merge the _SPECIES and _ALIAS sections of %registry_register.
    foreach my $section ( 'Species', 'Alias' ) {
      my $section_key = '_' . uc($section);

      while ( my ( $key, $value ) =
        each( %{ $registry_register{$section_key} } ) )
      {
        if ( !exists( $merged_register{$section_key}{$key} ) ) {
          $merged_register{$section_key}{$key} = $value;
        } elsif ($verbose) {
          printf( "%s '%s' found on multiple servers, "
              . "using first found\n",
            $section, $key );
        }
      }
    }
  } ## end foreach my $arg (@args)

  # Add the DBAs from the _SPECIES section into the _DBA section.
  foreach my $species_hash ( values( %{ $merged_register{_SPECIES} } ) )
  {
    foreach my $group_hash ( values( %{$species_hash} ) ) {
      if ( ref($group_hash) eq 'HASH' && exists( $group_hash->{_DB} ) )
      {
        push( @{ $merged_register{_DBA} }, $group_hash->{_DB} );
      }
    }
  }

  %registry_register = %merged_register;
  
  my $count = $self->get_DBAdaptor_count() - $original_count;
  return $count >= 0 ? $count : 0; 
} ## end sub load_registry_from_multiple_dbs

#
# Web specific routines
#

# =head2 DEPRECATED load_registry_with_web_adaptors
# 
#   DEPRECATED: Use load_registry_from_db instead.
# 
# =cut
# 
# sub load_registry_with_web_adaptors{
#   my $self = shift;
# 
#   deprecate('Use the load_registry_from_db instead'); 
#   my $site_eval = eval{ require SiteDefs };
#   if ($@ or (!defined($site_eval))){ die "Can't use SiteDefs.pm - $@\n"; }
#     SiteDefs->import(qw(:ALL));
# 
#   my $species_eval = eval{ require SpeciesDefs };
#   if ($@ or (!defined($species_eval))){ die "Can't use SpeciesDefs.pm - $@\n"; }
#   my $conf = new SpeciesDefs();
#   
#   my %species_alias = %{$SiteDefs::ENSEMBL_SPECIES_ALIASES};
# 
#   foreach my $spec (keys %species_alias){
#     Bio::EnsEMBL::Registry->add_alias($species_alias{$spec},$spec);
#   }
#   return;
# }
# 
# =head2 set_default_track
# 
#   Sets a flag to say that that this species/group are a default track and do not
#   need to be added as another web track.
# 
#   Arg [1]    : name of the species to get the adaptors for in the registry.
#   Arg [2]    : name of the type to get the adaptors for in the registry.
#   Example    : $merged = Bio::EnsEMBL::Registry->set_default_track("Human","core");
#   Returntype : none
#   Exceptions : none
#   Status     : At Risk.
# 
# =cut
# 
# sub set_default_track {
#   my ( $self, $species, $group ) = @_;
# 
#   $species = get_alias($species);
#   $self->storage->{'def_track'}{$species}{ lc($group) } = 1;
#   return;
# }
# 
# =head2 default_track
# 
#   Check flag to see if this is a default track
# 
#   Arg [1]    : name of the species to get the adaptors for in the registry.
#   Arg [2]    : name of the type to get the adaptors for in the registry.
#   Example    : $merged = Bio::EnsEMBL::Registry->set_default_track("Human","core");
#   Returntype : int 
#   Exceptions : none
#   Status     : At Risk.
# 
# =cut
# 
# sub default_track {
#   my ( $self, $species, $group ) = @_;
# 
#   $species = get_alias($species);
#   if (
#     defined( $self->storage->{'def_track'}{$species}{ lc($group) } ) )
#   {
#     return 1;
#   }
# 
#   return 0;
# }
# 
# 
# =head2 add_new_tracks
# 
#   Will add new gene tracks to the configuration of the WEB server if they are
#   not of the type default and the configuration already has genes in the display.
# 
#   Arg [1]    : hash of the default configuration of the web page
#   Returntype : none
#   Exceptions : none
#   Called by  : UserConfig.pm
#   Status     : At Risk.
#   
# =cut
# 
# sub add_new_tracks{
#   my($self, $conf, $pos) = @_;
# 
#   my $start = 0;
#   my $species_reg = $self->get_alias($conf->{'species'},"nothrow");
#   my %pars;
# #  print STDERR "Species $species_reg check for default tracks\n";
#   if(defined($species_reg)){
#     foreach my $dba (@{$self->get_all_DBAdaptors()}){
#       if(!$self->default_track($dba->species,$dba->group)){
#         $pars{'available'} = "species ".$self->get_alias($dba->species());
#         $pars{'db_alias'} = $dba->group();
# #       print STDERR "Adding new track for ".$dba->species."\t".$dba->group."\n";
#         $conf->add_new_track_generictranscript('',$dba->group(), "black",$pos,%pars);
#         $pos++;
#       }
#     }
#   }
#   return $pos;
# 
# }

=head2 no_version_check
  
  getter/setter for whether to run the version checking
  
  Arg[0]     : (optional) int
  Returntype : int or undef if not set
  Exceptions : none
  Status     : At Risk.

=cut
  
sub no_version_check {
  my ( $self, $arg ) = @_;
  ( defined $arg )
    && ( $self->storage->{'_no_version_check'} = $arg );

  return $self->storage->{'_no_version_check'};
}

=head2 no_cache_warnings

  Arg[0]      : boolean for turning the flag on and off
  Description : Turns off any warnings about not using caching in all available 
                adaptors. 
  Returntype  : boolean Current status
  Exceptions  : None

=cut

sub no_cache_warnings {
  my ($self, $arg) = @_;
  if(defined $arg) {
    $Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SILENCE_CACHE_WARNINGS = $arg;
  }
  return $Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SILENCE_CACHE_WARNINGS;
}

  
=head2 version_check
  
  run the database/API code version check for a DBAdaptor
  
  Arg[0]     : DBAdaptor to check
  Returntype : int 1 if okay, 0 if not the same 
  Exceptions : none
  Status     : At Risk.

=cut
  
  
sub version_check {
  my ( $self, $dba ) = @_;

  # Check the datbase and versions match
  # give warning if they do not.
  my $check = no_version_check();

  if ( (
      defined( $ENV{HOME} )
      and ( -e $ENV{HOME} . "/.ensemblapi_no_version_check" ) )
    or ( defined($check) and ( $check != 0 ) ) )
  {
    return 1;
  }

  my $mca =
    $self->get_adaptor( $dba->species(), $dba->group(),
    "MetaContainer" );

  my $database_version = 0;
  if ( defined($mca) ) {
    $database_version = $mca->get_schema_version();
  }

  if ( $database_version == 0 ) {
    # Try to work out the version
    if ( $dba->dbc()->dbname() =~ /^_test_db_/x ) {
      return 1;
    }
    if ( $dba->dbc()->dbname() =~ /(\d+)_\S+$/x ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_compara_(\d+)/x ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_help_(\d+)/x ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_ontology_(\d+)/x ) {
      $database_version = $1;
    } elsif ( $dba->dbc()->dbname() =~ /ensembl_stable_ids_(\d+)/x ) {
      $database_version = $1;
    } else {
      warn(
        sprintf(
          "No database version for database %s "
            . ". You must be using a post version 34 database "
            . "with version 34 or later code.\n"
            . "You need to update your database "
            . "or use the appropriate Ensembl software release "
            . "to ensure your script does not crash\n",
          $dba->dbc()->dbname() ) );
    }
  } ## end if ( $database_version...

  if ( $database_version != software_version() ) {
    warn(
      sprintf(
        "For %s there is a difference in the software release (%s) "
          . "and the database release (%s). "
          . "You should update one of these to ensure that your script "
          . "does not crash.\n",
        $dba->dbc()->dbname(),
        software_version(), $database_version
      ) );
    return 0;
  }

  return 1;    # Ok
} ## end sub version_check

=head2 get_all_species 

  Arg [1]    : String group type, such as core, or otherfeatures
  Description: Method for getting all valid species names found in available
               databases. This excludes the ancestral sequence databases, and
               any species from a non-core database. Specifying a group allows
               the list to apply to non-core database types.
  Example    : my @species_names = @{ $reg->get_all_species() };
  Returntype : Listref of species names
  
=cut

sub get_all_species {
    my ($self,$group) = @_;
    $group ||= 'core';
    my @species;
    my $storage = $self->storage();
    foreach my $name (keys %{$storage->{_SPECIES}}) {
        push @species, $name if (
            # limit species names to given db group and no ancestral dbs
            $storage->{_SPECIES}->{$name}->{$group}
            && $name !~ /^ancestral/i 
        );
    }
    return \@species;
}


=head2 get_species_and_object_type

  Description:  Get the species name, object type (gene, transcript,
                translation, or exon etc.), and database type for a
                stable ID.

  Arg [1]    :  String stable_id
                The stable ID to find species and object type for.

  Arg [2]    :  String known_type (optional)
                The type of the stable ID, if it is known.

  Arg [3]    :  String known_species (optional)
                The species, if known

  Arg [4]    :  String known_db_type (optional)
                The database type, if known
                
  Example    :  my $stable_id = 'ENST00000326632';

                my ( $species, $object_type, $db_type ) =
                  $registry->get_species_and_object_type($stable_id);

                my $adaptor =
                  $registry->get_adaptor( $species, $db_type,
                                          $object_type );

                my $object = $adaptor->fetch_by_stable_id($stable_id);

  Return type:  Array consisting of the species name, object type,
                and database type.  The array may be empty if no
                match is found.

  Exceptions :  none
  Status     :  At Risk.

=cut

sub get_species_and_object_type {
  my ($self, $stable_id, $known_type, $known_species, $known_group, $force_long_lookup) = @_;
  
  #Convert the known species into the real name if it was given
  $known_species = $self->get_alias($known_species) if $known_species;
  
  #Use the stable id if it was available & we did not force using long lookup
  my $stla = $self->get_adaptor('multi', 'stable_ids', 'stableidlookup');
  if($stla && ! $force_long_lookup) {
    return $stla->fetch_object_location($stable_id, $known_type, $known_species, $known_group);
  }
  
  #Scan using long lookup
  my $scanner = Bio::EnsEMBL::Utils::StableIdsDBAdaptorLookup->new($self);
  return $scanner->fetch_object_location($stable_id, $known_type, $known_species, $known_group);
}

1;