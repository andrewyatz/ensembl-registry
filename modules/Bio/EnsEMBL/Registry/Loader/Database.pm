package Bio::EnsEMBL::Registry::Loader::Database;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Registry::Loader::Base/;

use Bio::EnsEMBL::ApiVersion qw(software_version);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning);

my $PRODUCTION_NAME_META_KEY = 'species.production_name';
my $ALIAS_META_KEY = 'species.alias';

=head2 load_registry

  Arg [HOST] : string
                The domain name of the database host to connect to.

  Arg [PORT] : (optional) integer
                The port to use when connecting to the database.

  Arg [USER] : string
                The name of the database user to connect with.

  Arg [PASS] : (optional) string
                The password to be used to connect to the database.

  Arg [VERBOSE]: (optional) boolean
                Whether to print database messages. This includes a listing
                of all available species & databases.

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
  my ( $host,         $port,     $user,
       $pass,         $verbose,  $db_version,
       $wait_timeout, $no_cache, $species, $species_suffix )
    = rearrange( [ 'HOST',          'PORT',
                   'USER',          'PASS',
                   'VERBOSE',       'DB_VERSION',
                   'WAIT_TIMEOUT',  'NO_CACHE',
                   'SPECIES',       'SPECIES_SUFFIX' ],
                 @args );

  
  my $registry = $self->registry();
  
  #If we were given a -VERBOSE or -NO_CACHE flag then set them & grab
  #back since it could have already been set
  $self->verbose($verbose) if $verbose;
  $verbose = $self->verbose($verbose);
  $self->no_cache($no_cache) if $no_cache;
  $no_cache = $self->no_cache();

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
    -NO_CACHE     => $no_cache
  );
  
  #Single DBC for all DB operations. Force a connection by requesting the handle
  #and causing an early bail from here rather than lower down the stack
  my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%basic_args);
  $dbc->db_handle();
  
  my $regex_patterns = $self->_patterns();
  my $alias_available = $self->_alias_available();
  my $default_aliases = $self->_default_aliases($species_suffix);
  my $species_filter_groups = $self->_species_filter_groups();
  
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
      
      my $parsed_species = undef;
      my $database_version = undef;
      my $multispecies = 0;
      
      # If we have a collection regex & it matched then we have a multi-species
      # database.
      if(defined $collection_re && $db_name =~ $collection_re) {
        ($parsed_species, $database_version) = ($1, $2);
        $multispecies = 1;
      }
      #Otherwise just a normal DB
      if(! $multispecies && $db_name =~ $single_re) {
        ($parsed_species, $database_version) = ($1, $2);
      }
      
      #Skip if we had no hit for species. Means we do not want to process this for this group
      next if ! defined $parsed_species;
      
      #Check if it matches the $parsed_species matches the $species value
      next if defined $species && $species_filter_groups->{$group} && $parsed_species !~ /^$species/;
      
      #Process only if the DB version matched the requested
      if($version == $database_version) {
        #Convert from a name & a connection to a 2D array of [species_id,name]
        my $to_load = $self->_get_dbas_to_load($dbc, $db_name, $parsed_species, $multispecies);
        my @created_dbas;
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
          push(@created_dbas, $dba);
          
          printf( "Species '%s' loaded from database '%s' (group : %s | id : %d)\n", $target_species, $db_name, $target_group, $target_species_id ) if $verbose;
        }

        #Apologies this is nasty. If we have a DB type that supports aliases then attempt a load.
        if($alias_available->{$group}) {
          # If it was a multi-species DB then do a batch load for all aliases (for speeds sake)
          if($multispecies) {
            my $aliases = $self->get_all_database_aliases($dbc, $db_name, $species_suffix);
            foreach my $parsed_species (keys %{$aliases}) {
              $registry->add_alias($parsed_species, $aliases->{$species});
            }
          }
          #Otherwise we just do one load for the single DB (there should only be one dba created)
          else {
            my $dba = shift @created_dbas;
            my $aliases = $self->get_aliases($dbc, $dba, $species_suffix);
            $registry->add_alias($dba->species(), $aliases);
          }
        }
      }
      else {
        printf("Not processing '%s' (%s) as its detected version %d is not the same as the requested version %d\n", 
          $parsed_species, $db_name, $database_version, $version) if $verbose;        
      }
      
      splice(@databases, $i, 1); #remove current element
      if($#databases == -1) { # if the last element index is -1 then finish the loop; the array is exhausted
        last;
      }
      #otherwise bring the iterator back one and decrement the size of the array by one
      $i--;
      $db_count--;
    }
    
    foreach my $key (keys %{$default_aliases}) {
      $registry->add_alias($key, $default_aliases->{$key});
    }
    
  }
  
  $dbc->disconnect_if_idle();
  return;
}

=head2 get_aliases

  Arg[1]     : (optional) Bio::EnsEMBL::DBSQL::DBConnection connection to the DB server. If
               not given we will use the supplied DBAdaptor's DBConnection
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
  $dbc ||= $dba->dbc();
  my ($quoted_table_name) = @{$dbc->quote_identifier([undef, $dba->dbc()->dbname(), 'meta'])};
  my $sql = sprintf('SELECT meta_value FROM %s WHERE meta_key =?', $quoted_table_name);
  my @params = ($ALIAS_META_KEY);
  if($dba->is_multispecies()) {
    $sql .= ' and species_id =?';
    push(@params, $dba->species_id());
  }
  my $meta_results = $dbc->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => \@params);
  return [ map { $_.$species_suffix } @{$meta_results}];
}

=head2 get_all_database_aliases

  Arg[1]     : Bio::EnsEMBL::DBSQL::DBConnection connection to the DB server
  Arg[2]     : (optional) String Database to look in. If not given then we use the database name
               the DBConnection instance is pointing to
  Arg[3]     : (optional) String Species suffix to add to all found aliases 
  Example    : $database_loader->get_all_database_aliases($dbc);
  Description: Queries the meta table for the given DBConnection. Should only be run
               on DBConnections which can provide this kind of table and data 
               (such as core). The database should provide the species name in the
               species.production_name meta value and aliases in the species.alias meta value
  Returntype : HashRef keyed by species name and value is an ArrayRef of aliases
  Exceptions : Thrown in the event of a DBI error
  Caller     : General
  Status     : At risk

=cut

sub get_all_database_aliases {
  my ($self, $dbc, $dbname, $species_suffix) = @_;
  $species_suffix ||= q{};
  $dbname ||= $dbc->dbname();
  my ($quoted_table_name) = @{$dbc->quote_identifier([undef, $dbname, 'meta'])};
  my $sql_template = <<'SQL';
SELECT m1.meta_value, m2.meta_value 
FROM %s m1 
JOIN %s m2 on (m1.species_id = m2.species_id) 
WHERE m1.meta_key =? AND m2.meta_key =?
SQL
  my $sql = sprintf($sql_template, $quoted_table_name, $quoted_table_name);
  my $params = [$PRODUCTION_NAME_META_KEY, $ALIAS_META_KEY];
  my $meta_results = $dbc->sql_helper()->execute_into_hash(-SQL => $sql, -PARAMS => $params, -CALLBACK => sub {
    my ( $row, $value ) = @_;
    if ( defined $value ) {
      push( @{$value}, $row->[1] );
      return;
    }
    my $new_value = [ $row->[1] ];
    return $new_value;
  });
  return $meta_results;
}

=head2 _get_version

  Arg[1]      : (optional) Integer The possible db_version
  Description : Provides a default version which is software version
                if a db_version was not defined
  Returntype  : Integer The version of Ensembl to query for

=cut

sub _get_version {
  my ($self, $db_version) = @_;
  return (defined $db_version) ? $db_version : software_version();
}

=head2 _get_databases

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection Used to query a server with
  Arg[2]      : Integer The version of Ensembl to look for
  Description : Queries the provided database instance for all databases matching
                the given version or the userdata DBs. All DBs should be
                post-processed to ensure they match the correct release version
  Returntype  : ArrayRef of possible database names

=cut

sub _get_databases {
  my ($self, $dbc, $version) = @_;
  my @dbs;
  foreach my $pattern (("%\\_${version}%",'userdata%')) {
    push(@dbs, @{$dbc->sql_helper()->execute_simple(-SQL => 'SHOW DATABASES LIKE ?', -PARAMS => [$pattern])});
  }
  return \@dbs;
}

=head2 _get_dbas_to_load

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection Used to query a server with
  Arg[2]      : String The database name to query against
  Arg[3]      : String Potential species name
  Arg[4]      : Boolean Flags if this is a multi-species database
  Description : Used to direct queries off to C<_multispecies_species()> or to 
                return the species with a species_id of 1.
  Returntype  : ArrayRef (2D). Structure is [[species_id, species_name]]

=cut

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

=head2 _multispecies_species

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection Used to query the server with
  Arg[2]      : String The database name to query against
  Description : Query a schema for all available multi-species which is done by 
                querying for all production names and their associcated species identifier
  Returntype  : ArrayRef (2D). Structure is [[species_id, species_name]]

=cut

sub _multispecies_species {
  my ($self, $dbc, $db_name) = @_;
  my ($table) = @{$dbc->quote_identifier([undef, $db_name, 'meta'])};
  my $sql = sprintf('SELECT species_id, meta_value FROM %s WHERE meta_key =?', $table);
  return $dbc->sql_helper()->execute(-SQL => $sql, -PARAMS => [$PRODUCTION_NAME_META_KEY]);
}

=head2 _species_filter_groups

  Description : Defines a lookup Hash of groups where you can apply 
                species filtering to respect the -SPECIES flag
  Returntype  : HashRef key is the group and value is a boolean 

=cut

sub _species_filter_groups {
  my ($self) = @_;
  return {
    core => 1,
    otherfeatures => 1,
    cdna => 1,
    vega => 1,
    rnaseq => 1,
    variation => 1,
    funcgen => 1,
    userupload => 1,
  };
}

=head2 _database_order

  Description : Defines the order in which we process groups. They are
                core, core-like, variation, funcgen, userupload,
                compara, ancestral, ontology and stable_ids
  Returntype  : ArrayRef containing Strings of the aformentioned groups
=cut

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

=head2 _alias_available

  Description : The database groups aliases can be retrieved from.
                Defaults to core and compara
  Returntype  : HashRef of groups. Keys are the group names

=cut

sub _alias_available {
  my ($self) = @_;
  return {
    core => 1,
    compara => 1
  };
}

=head2 _default_aliases

  Arg[1]      : (optional) String Apply a species suffix to all names and aliases
  Description : Provides a number of default aliases. Notably it adds
                multi -> compara, ontology and stable_ids
                Ancestral sequences -> ancestral_sequences
  Returntype  : HashRef keyed by the species name. Value is an ArrayRef of aliases

=cut

sub _default_aliases {
  my ($self, $species_suffix) = @_;
  $species_suffix ||= q{};
  return {
    "multi${species_suffix}" => [
      "compara${species_suffix}", "ontology${species_suffix}", "stable_ids${species_suffix}"
    ],
    "Ancestral sequences${species_suffix}"  => [
      "ancestral_sequences${species_suffix}"
    ],
  };
}

=head2 _post_process_name

  Arg[1]      : String The group to process
  Arg[2]      : String The species name
  Description : Provides some ability to rename a database species 
                because its offical registry name differs 
                from that which was picked up from the DB name
  Returntype  : String The processed name to use
=cut

sub _post_process_name {
  my ($self, $group, $species_name) = @_;
  my $altered_name = $species_name;
  if($group eq 'compara') {
    #Regex which captures any additional names given to a DB e.g. pan_homology or plants
    if($species_name =~ $self->_patterns()->{NAMES}->{compara}) {  
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

=head2 _post_process_group

  Arg[1]      : String The group to process
  Description : Force group changes if needed. Ancestral (detected group) should 
                be core once in the registry. The group ancestral comes from an
                internal need to know when to grep with a different regex
  Returntype  : String Processed group to use

=cut

sub _post_process_group {
  my ($self, $group) = @_;
  if($group eq 'ancestral') {
    return 'core';
  }
  return $group;
}

=head2 _patterns

  Description : Provides a number of pre-compiled regular expressions. 
                Cached on $self for convenience and to avoid overheads 
                of repeated calling. This should not be a huge issue IMHO
  Returntype  : HashRef of all available regex patterns

=cut

sub _patterns {
  my ($self) = @_;
  
  return $self->{regex_patterns} if exists $self->{regex_patterns};
    
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
    },
    
    END         =>  qr/(?:_\d+)? _ (\d+) _ \d+? \w? /xms, # Represents [_10]_(63)_1[a] where the E! version is captured. The \d+? is important to capture assemblies like 235
  );
    
  my $regex_patterns = {
    PARTIAL => \%RE,    # Provide subcallers with the partial RegEx to avoid re-implementation
    NAMES => {          # Give the ability to search for a name and capture it
      compara => qr/^ensembl_compara_ (\w+?) (?:_\d+)? _ \d+ $/xms, #tested with ensembl_compara_72, ensembl_compara_plants_10_72 and ensembl_compara_pan_homology_10_72
    },
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
  
  return $self->{regex_patterns} = $regex_patterns;
}

1;