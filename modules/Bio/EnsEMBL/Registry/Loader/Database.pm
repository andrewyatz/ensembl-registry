package Bio::EnsEMBL::Registry::Loader::Database;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Registry::Loader::Base/;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning);

=head2 load_registry

  Arg [HOST] : string
                The domain name of the database host to connect to.

  Arg [USER] : string
                The name of the database user to connect with.

  Arg [PASS] : (optional) string
                The password to be used to connect to the database.

  Arg [PORT] : (optional) integer
                The port to use when connecting to the database.

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

   Arg [SPECIES_SUFFIX]: (optional) string
                This option will append the string to the species name
                in the registry for all databases found on this server.

  Example :

    $database_loader->load_registry(
      -host    => 'ensembldb.ensembl.org',
      -user    => 'anonymous',
      -verbose => '1'
    );

  Description: Will load the correct versions of the Ensembl
               databases for the software release it can find on a
               database instance into the registry.  Also adds a set
               of standard aliases.
  Returntype : Nothing
  Exceptions : Thrown if the given MySQL database cannot be connected to
               or there is any error whilst querying the database.
  Status     : Stable

=cut

sub load_registry {
  my ($self, @args) = @_;
  my (  $host,         $port,     $user,
        $pass,         $db_version,
        $wait_timeout, $species,
        $species_suffix )
    = rearrange( [  'HOST',          'PORT',
                    'USER',          'PASS',
                    'DB_VERSION',
                    'WAIT_TIMEOUT',
                    'SPECIES',       'SPECIES_SUFFIX' ], @args );

  my $registry = $self->registry();
  my $verbose = $self->verbose();
  my $no_caching = $self->no_caching();

  if(defined $species) {
    $species = lc($species);
    $species =~ tr/ -/__/;
  }
  $species_suffix = q{} if ! defined $species_suffix;
  
  if(! defined $db_version) {
    # Do checking for the -DB_VERSION flag which can be mis-spelt. Regex assembled using:
    # perl -MRegexp::Assemble -e '$r=Regexp::Assemble->new(); $r->add($_) for ("-dbversion","-version","-verion","-verison"); print $r->re, "\n";'
    my %hashed_args = @args;
    my ($possible_key) = grep { $_ =~ /(?-xism:-(?:ver(?:is?|si)|dbversi)on)/xism } keys %hashed_args;
    if($possible_key) {
      my $msg = sprintf(q{Detected no -DB_VERSION flag but found '%s'; assuming a mis-spelling. Please fix}, $possible_key);
      warning($msg);
      $db_version = $hashed_args{$possible_key};
    }
  }
  
  $user ||= "ensro";
  if ( !defined($port) ) {
    $port = 3306;
    if ( $host eq "ensembldb.ensembl.org" && defined($db_version) && $db_version < 48 ) {
      $port = 4306;
    }
  }

  $wait_timeout ||= 0;
  my $version = (defined $db_version) ? $db_version : software_version();
  printf( "Will only load v%d databases\n", $version ) if $verbose;
  
  my %basic_args = (
    -HOST         => $host,
    -USER         => $user,
    -PASS         => $pass,
    -PORT         => $port,
    -WAIT_TIMEOUT => $wait_timeout,
    -NO_CACHE     => $no_caching
  );
  
  #Single DBC for all DB operations
  my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%basic_args);
  
  my $regex_patterns = $self->_patterns();
  my $alias_available = $self->_alias_available();
  my $default_aliases = $self->_default_aliases($species_suffix);
  
  #Query database for the available databases
  my @databases = @{$self->_get_databases($dbc, $version)};
  
  #Iterate through the database groups we need to load. 
  #core goes first followed by core-like, variation, regulation, userdata, compara & supplementals
  foreach my $group (@{$self->_database_order()}) {
    
    printf("Working with group '%s'\n", $group) if $verbose;
    
    #Attempt to import the group's default module. If this fails skip this database group
    my $module = $self->group_to_adaptor($group);
    my $success = $self->import_package($module, 0);

    #Failed so skip this group
    if(! $success) {
      printf("%s module not found on PERL5LIB so %s databases will be ignored\n", $module, $group) if $verbose;
      next;
    }
    
    my $regular_expressions = $regex_patterns->{FULL}->{$group};
    my $single_re = $regular_expressions->{single};
    my $collection_re = $regular_expressions->{collection};
    
    #Scan the entire database array per group; cannot assume any ordering
    my $db_count = scalar(@databases);
    for(my $i = 0; $i < $db_count; $i++) {
      my $db_name = $databases[$i];
      
      #happens because we undef the DB names once we've successfully
      #processed it. Sucky logic.
      #TODO come up with some better logic to process a list like a Java
      #iterator e.g. iter.remove() the current element from a list & then call iter.next()
      next unless defined $db_name;
      
      my $species = undef;
      my $database_version = undef;
      my $multispecies = 0;
      
      # If we have a collection regex & it matched then we have a multi-species
      # database.
      if(defined $collection_re && $db_name =~ $collection_re) {
        ($species, $database_version) = ($1, $2);
        $multispecies = 1;
      }
      #Otherwise just a normal DB
      if(! $multispecies && $db_name =~ $single_re) {
        ($species, $database_version) = ($1, $2);
      }
      
      #Skip if we had no hit for species. Means we do not want to process this for this group
      next if ! defined $species;
      
      #Process only if the DB version matched the requested
      if($version == $database_version) {
        #Convert from a name & a connection to a 2D array of [species_id,name]
        my $to_load = $self->_get_dbas_to_load($dbc, $db_name, $species, $multispecies);
        foreach my $dba_to_load (@{$to_load}) {
          my ($target_species_id, $target_species) = @{$dba_to_load};
          
          #Make sure we've got the right name and group
          $target_species = $self->_post_process_name($group, $target_species);
          my $target_group = $self->_post_process_group($group);
          
          #Create the DBAdaptor, find the aliases and add them in
          my $dba = $module->new(
            -SPECIES => $target_species.$species_suffix,
            -SPECIES_ID => $target_species_id,
            -GROUP => $target_group,
            -IS_MULTISPECIES => $multispecies,
            -DBNAME => $db_name,
            %basic_args
          );
          
          my $aliases = ($alias_available->{$group}) ? $self->get_aliases($dbc, $dba, $species_suffix) : [] ; 
          $registry->add_alias($species, $aliases);
          
          printf( "Species '%s' (id:%d) group '%s' loaded from database '%s'\n", $target_species, $target_species_id, $target_group, $db_name ) if $verbose;
        }
      }
      else {
        printf("Not processing '%s' (%s) as its detected version %d is not the same as the requested version %d\n", 
          $species, $db_name, $database_version, $version) if $verbose;        
      }
      
      $databases[$i] = undef;
    }
    
    foreach my $key (keys %{$default_aliases}) {
      $registry->add_alias($key, $default_aliases->{$key});
    }
    
  }
  
  $dbc->disconnect_if_idle();
  return;
}

=head2 get_aliases

  Arg[1]     : Bio::EnsEMBL::DBSQL::DBConnection connection to the DB server
  Arg[2]     : Bio::EnsEMBL::DBSQL::DBAdaptor instance/species to lookup aliases for
  Arg[3]     : (optional) String Species suffix to add to all found aliases 
  Example    : $database_loader->get_aliases($dbc, $dba);
  Description: Queries the meta table for the given DBAdaptor (using the provided
               DBConnection) for all available aliases. Should only be run
               on DBAdaptors which can provide this kind of table and data 
               (such as core and compara)
  Returntype : ArrayRef of all known aliases for this species
  Exceptions : Thrown in the event of a DBI error
  Caller     : General
  Status     : At risk

=cut

sub get_aliases {
  my ($self, $dbc, $dba, $species_suffix) = @_;
  $species_suffix ||= q{};
  my ($quoted_table_name) = $dbc->quote_identifier([undef, $dba->dbc()->dbname(), 'meta']);
  my $sql = sprintf('SELECT meta_value FROM %s WHERE meta_key =? and species_id =?', $quoted_table_name);
  my $meta_results = $dbc->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => ['species.alias', $dba->species_id()]);
  return [ map { $_.$species_suffix } @{$meta_results}];
}

sub _get_version {
  my ($self, $db_version) = @_;
  return (defined $db_version) ? $db_version : software_version();
}

sub _get_databases {
  my ($self, $dbc, $version) = @_;
  my @dbs;
  foreach my $pattern (("%\\_${version}%",'userdata%')) {
    push(@dbs, @{$dbc->sql_helper()->execute_simple(-SQL => 'SHOW DATABASES LIKE ?', -PARAMS => [$pattern])});
  }
  return \@dbs;
}

sub _get_dbas_to_load {
  my ($self, $dbc, $db_name, $species, $multispecies) = @_;
  my $output;
  if($multispecies) {
    $output = $self->_multispecies_species($dbc, $db_name);
  }
  else {
    $output = [[1, $species]];
  }
  return $output;
}

sub _multispecies_species {
  my ($self, $dbc, $db_name) = @_;
  my ($table) = $dbc->quote_identifier([undef, $db_name, 'meta']);
  my $sql = sprintf('SELECT species_id, meta_value FROM %s WHERE meta_key =?', $table);
  return $dbc->sql_helper()->execute(-SQL => $sql, -PARAMS => ['species.db_name']);
}

sub _database_order {
  my ($self) = @_;
  return [qw/
    core otherfeatures cdna vega rnaseq 
    variation 
    funcgen
    userupload
    compara ancestral
    ontology
    stable_ids
  /];
}

sub _alias_available {
  my ($self) = @_;
  return {
    map { $_ => 1 }
    qw/core compara/
  };
}

sub _default_aliases {
  my ($self, $species_suffix) = @_;
  $species_suffix ||= q{};
  return {
    "multi${species_suffix}" => [
      "compara${species_suffix}", "ontology${species_suffix}", "stable_ids${species_suffix}"
    ],
    "Ancestral sequences${species_suffix}"  => [
      "ancestal_sequences${species_suffix}"
    ],
  };
}

#Provides some ability to rename a database species because its offical
#registry name differs from that which was picked up from the DB name
sub _post_process_name {
  my ($self, $group, $species_name) = @_;
  my $altered_name = $species_name;
  if($group eq 'compara') {
    #inline regex required due to different capture requirements from the defaults
    if($species_name =~ /^ensembl_compara_([a-zA-Z_])_\d+/xms) { 
      $altered_name = $1;
    }
    else {
      $altered_name = 'multi';
    }
  }
  elsif($group eq 'ontology' || $group eq 'stable_ids') {
    $altered_name = 'multi';
  }
  elsif($group eq 'ancestral') {
    $altered_name = 'Ancestral sequences';
  }
  return $altered_name;
}

#Force group changes if needed. Ancestral (detected group) should 
#be core once in the registry. The group ancestral comes from an
#internal need to know when to grep with a different regex
sub _post_process_group {
  my ($self, $group) = @_;
  if($group eq 'ancestral') {
    return 'core';
  }
  return $group;
}

sub _patterns {
  my ($self) = @_;
    
  my %RE = (
    NAME        => qr{ [a-z]+ _ [a-z0-9]+ (?:_[a-z0-9]+)? }xms, #Name of the DB (binomial & trinomial accepted)
    
    COLLECTION  => qr{ \w+_collection }xms,     #Indicates multispecies DB
    
    E_VERSION   => qr{ \d+ },                   #EnsemblGenomes optional version
    EG_VERSION  => qr{ (?:_\d+)? },             #Ensembl version
    
    DB_VERSION  => qr{ _\w+ }xms,               #End of the name
    
    GROUPS => {                                 #All known groups
      CORE            => qr{_core           }xms,
      CDNA            => qr{_cdna           }xms,
      VEGA            => qr{_vega           }xms,
      OTHER_FEATURES  => qr{_otherfeatures  }xms,
      RNASEQ          => qr{_rnaseq         }xms,
      USER_UPLOAD     => qr{_userdata       }xms,
      VARIATION       => qr{_variation      }xms,
      FUNCGEN         => qr{_funcgen        }xms,
      COMPARA         => qr{_compara        }xms,
      ANCESTRAL       => qr{_ancestral      }xms,
      ONTOLOGY        => qr{_ontology       }xms,
      STABLE_IDS      => qr{_stable_ids     }xms,
    }
  );
  
  # Represents [_10]_(63)_1[a] where the E! version is captured
  $RE{END} = qr/(?:_\d+)? _ (\d+) _ \d \w? /xms;
    
  return {
    PARTIAL => \%RE,    # Provide subcallers with the partial RegEx to avoid re-implementation
    FULL    => {        # All full regular expressions capture a name and a E! version
      core => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{CORE} $RE{END} $/xms,       # homo_sapiens_core_63_37
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{CORE} $RE{END} $/xms    # eschericia_shigella_collection_core_10_63_1
      },
      cdna => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{CDNA} $RE{END} $/xms,       # homo_sapiens_cdna_63_37
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{CDNA} $RE{END} $/xms    # eschericia_shigella_collection_cdna_10_63_1
      },
      vega => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{VEGA} $RE{END} $/xms,
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{VEGA} $RE{END} $/xms
      },
      otherfeatures => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{OTHER_FEATURES} $RE{END} $/xms,
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{OTHER_FEATURES} $RE{END} $/xms
      },
      rnaseq => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{RNASEQ} $RE{END} $/xms,
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{RNASEQ} $RE{END} $/xms
      },
      variation => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{VARIATION} $RE{END} $/xms,
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{VARIATION} $RE{END} $/xms
      },
      funcgen => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{FUNCGEN} $RE{END} $/xms,
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{FUNCGEN} $RE{END} $/xms
      },
      userupload => {
        single      => qr/^ ($RE{NAME})   $RE{GROUPS}{USER_UPLOAD} $/xms,             #homo_sapiens_userdata
        collection  => qr/^ ($RE{COLLECTION}) $RE{GROUPS}{USER_UPLOAD} $/xms          #escherichia_shigella_collection_userdata
      },
      compara => {
        single      => qr/^ (ensembl_compara (?:_\w+)?) (?:_\d+)? _ (\d+) $/xms,      #ensembl_compara_63 and ensembl_compara_fungi_10_63
      },
      ancestral => {
        single      => qr/^ (ensembl_ancestral) (?: _ \w+)? (?:_\d+)? _ (\d+) $/xms,  #ensembl_ancestral_63 and ensembl_ancestral_eg_10_63
      },
      ontology => {
        single      => qr/^ (ensembl(?:[a-z]+)?_ontology) (?:_\d+)? _ (\d+) $/xms,    #ensembl_ontology_63 and ensemblgenomes_ontology_10_63
      },
      stable_ids => {
        single      => qr/^ (ensembl(?:[a-z]+)?_stable_ids) (?:_\d+)? _ (\d+) $/xms,  #ensembl_stable_ids_63 and ensemblgenomes_stable_ids_10_63
      },
    }
  };
}

1;