package Bio::EnsEMBL::Registry::Loader::PerlHash;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Registry::Loader::Base/;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

=head2 load_registry

  Arg [1]    : HashRef $perl_hash raw HASH reference to load a Registry from. Hash requires
               two keys. First is adaptors pointing at an array of hash refs. Keys of these
               hash refs should be the same as construction parameters to DBAdaptor minus
               any rearrange semantics and lowercased. Second key is aliases and should
               point to a hash keyed by the known species name and then an array of aliases.
               See the example for a concrete example of the format.
  Example    : my $hash = { 
                adaptors => [
                  { host => 'host', port=>3306, user=>'anon', pass=>'pass', dbname=>'db', species => 'homo_sapiens' },
                  { host => 'host', port=>3306, user=>'anon', pass=>'pass', dbname=>'another_db', species => 'mouse' },
                ], 
                aliases => { 'homo_sapiens' => ['human'], 'mouse' => ['squeek']}
               };
               $perl_hash_loader->load_registry($hash);
  Description: Parses the given HashRef and loads the
               associcated Registry with the specified adaptors
  Returntype : None
  Exceptions : Thrown if $hash was not defined or not a HashRef
  Caller     : General
  Status     : At risk

=cut

sub load_registry {
  my ($self, $perl_hash) = @_;
  assert_ref($perl_hash, 'HASH', 'perl_hash');
  my $registry = $self->registry();
  my $verbose = $self->verbose();
  my $no_cache = $self->no_cache();
  
  my %skip_group;
  if(exists $perl_hash->{adaptors}) {
    my $adaptors = $perl_hash->{adaptors};
    assert_ref($adaptors, 'ARRAY', 'adaptors');
    my $adaptor_count = 1;
    foreach my $adaptor_hash (@{$adaptors}) {
      #Find the group, lookup its DBAdaptor and import. 
      #Skip all DBAdaptors linked to this group if any stage fails
      my $group = $adaptor_hash->{group};
      if(! $group) {
        warning sprintf("Key 'group' is missing for adaptor specification %d. Skipping.\n", $adaptor_count);
        next;
      }
      next if $skip_group{$group};
      my $module = $self->group_to_adaptor($group);
      
      if (!defined($module)) {
        $skip_group{$group} = 1;
        warning sprintf("Group '%s' does not convert to a known module. Skipping\n", $group);
        next;
      }
      
      #Import the module for this group if we need to
      if(!$self->import_package($module)) {
        $skip_group{$group} = 1;
        #Not an error; if the user has not imported Variation, Compara or Regulation then the
        #import will legitimatley fail
        printf("Cannot load group '%s' as the module '%s' cannot be imported.\n", $group, $module) if $verbose;
        next;
      }
      
      printf("Loading group '%s' with module '%s'\n", $group, $module) if $verbose;
      
      #Convert from perl hash into one more compatible with DBAdaptor and rearrange
      my %adaptor_args;
      foreach my $parameter ( keys %{$adaptor_hash} ) {
        $adaptor_args{"-${parameter}"} = $adaptor_hash->{$parameter};
      }
      
      my $dba = $module->new(%adaptor_args, -NO_CACHE => $no_cache, -REGISTRY => $registry);
      
      #Report the load
      printf( "Species '%s' (id:%d) group '%s' loaded\n", $dba->species(), $dba->species_id(), $group) if $verbose;
      
      $adaptor_count++;
    }
    print("Loaded species\n") if $verbose;
  }
  
  if(exists $perl_hash->{aliases}) {
    my $aliases_hash = $perl_hash->{aliases};
    assert_ref($aliases_hash, 'HASH', 'aliases');
    foreach my $species (keys %{$aliases_hash}) {
      my $aliases = $perl_hash->{aliases}->{$species};
      $registry->add_alias($species, $aliases)
    }
    print("Loaded aliases\n") if $verbose;
  }
  return;
}

=head2 serialise_registry

  Example    : my $hash = $perl_hash_loader->serialise_registry();
  Description: For the associated registry object this methods returns the HashRef
               representation
  Returntype : HashRef of a serialised registry
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub serialise_registry {
  my ($self) = @_;
  my $registry = $self->registry();
  
  my $hash = {
    adaptors => [],
    aliases => {},
  };
  
  my %species_names;
  
  my @dbas = @{$registry->get_all_DBAdaptors()};
  foreach my $dba (@dbas) {
    $species_names{$dba->species()} = 1;
    my $group = $dba->group();
    my $dbc = $dba->dbc();
    
    my $adaptor = {
      species => $dba->species(),
      group => $group,
      host => $dbc->host(),
      port => $dbc->port(),
      user => $dbc->username(),
      pass => $dbc->password(),
      dbname => $dbc->dbname(),
      driver => $dbc->driver(),
    };
    
    $adaptor->{disconnect_when_inactive} = $dbc->disconnect_when_inactive() if $dbc->disconnect_when_inactive();
    $adaptor->{wait_timeout} = $dbc->timeout() if $dbc->timeout();
    $adaptor->{reconnect_when_lost} = $dbc->reconnect_when_lost() if $dbc->reconnect_when_lost();
    
    if($dba->is_multispecies()) {
      $adaptor->{multispecies_db} = 1;
      $adaptor->{species_id} = $dba->species_id();
    }
    
    push(@{$hash->{adaptors}}, $adaptor);
  }
  
  foreach my $species (keys %species_names) {
    my $aliases = $registry->get_all_aliases($species);
    next unless @{$aliases};
    #sort aliases into alphabetical order. Cheap win since C sort (no closure) is very fast
    $hash->{aliases}->{$species} = [sort @{$aliases}]; 
  }
  
  return $hash;
}

1;