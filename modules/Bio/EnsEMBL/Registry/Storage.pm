package Bio::EnsEMBL::Registry::Storage;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref wrap_array);

##### Constructor

sub new {
  my ($class) = @_;
  $class ||= ref($class);
  my $self = bless({
    alias       => {},
    db_adaptor  => {},
    dna_adaptor => {}
  }, $class);
  return $self;
}

##### Retrieval methods

sub get_DBAdaptor {
  my ($self, @args) = @_;
  my ($species, $group, $no_warning) = rearrange([qw(species group no_warning)], @args);
  
  throw "Species was not defined" if ! defined $species;
  my $alias = $self->get_alias($species);
  
  my $dba = $self->{db_adaptor}->{$species}->{$group}; 
  
  if(! $no_warning) {
    warning "The species '$species' and group '$group' resulted in no found DBAdaptor";
  }
  
  return $dba;
}

sub get_all_DBAdaptors {
  my ($self, @args) = @_;
  my ($species, $group) = rearrange([qw(species group)], @args);
  
  my $adaptors;
  
  if(defined $species && exists $self->{db_adaptor}->{$species}) {
    $species = $self->get_alias($species);
    
    if(defined $group) {
      $adaptors = [$self->get_DBAdaptor(-SPECIES => $species, -GROUP => $group, -NO_WARNING => 1)];
    }
    else {
      $adaptors = [@{$self->{db_adaptor}->{$species}}];
    }
  }
  elsif(defined $group) {
    $adaptors = $self->grep_DBAdaptors(sub {
      my ($dba) = @_;
      return ( $dba->group() eq $group ) ? 1 : 0;
    });
  }
  else {
    $adaptors = $self->grep_DBAdaptors(sub {return 1;});
  }
  
  return $adaptors;
}

sub get_DNAAdaptor {
  my ($self, @args) = @_;
  my ($adaptor, $species, $group) = $self->_process_args([], @args);
  my $dna_args = $self->{dna_adaptor}->{$species}->{$group};
  return unless defined $dna_args;
  my ($dna_species, $dna_group) = @{$dna_args};
  return $self->get_DBAdaptor(-SPECIES => $dna_species, -GROUP => $dna_group);
}

sub get_alias {
  my ($self, $name) = @_;
  return $self->{alias}->{lc($name)};
}

sub get_all_aliases {
  my ($self, $name) = @_;
  my @aliases;
  my $species = $self->get_alias($name);
  foreach my $key ( %{$self->{alias}} ) {
    next if $species eq $key;
    my $target_species = $self->{alias}->{$key};
    push(@aliases, $key) if $species eq $target_species;
  }
  return \@aliases;
}

##### Query methods

sub alias_exists {
  my ($self, $name) = @_;
  return (exists $self->{alias}->{lc($name)}) ? 1 : 0;
}

sub DBAdaptor_exists {
  my ($self, @args) = @_;
  my ($adaptor, $species, $group) = $self->_process_args([], @args);
  return (exists $self->{db_adaptor}->{$species}->{$group}) ? 1 : 0;
}

sub grep_apply_to_DBAdaptors {
  my ($self, $filter, $apply) = @_;
  assert_ref($apply, 'CODE');
  my $adaptors = $self->grep_DBAdaptors($filter);
  foreach my $adaptor (@{$adaptors}) {
    $apply->($adaptor);
  } 
  return;
}

sub grep_DBAdaptors {
  my ($self, $callback) = @_;
  assert_ref($callback, 'CODE');
  my @adaptors;
  foreach my $species (keys %{$self->{db_adaptor}}) {
    foreach my $group (keys %{$self->{db_adaptor}->{$species}}) {
      my $adaptor = $self->{db_adaptor}->{$species}->{$group};
      push(@adaptors, $adaptor) if $callback->($adaptor);
    }
  }
  return \@adaptors;
}

##### Addition methods

sub add_DBAdaptor {
  my ($self, @args) = @_;
  my ($adaptor, $species, $group) = $self->_process_args([], @args);  
  assert_ref($adaptor, 'Bio::EnsEMBL::DBSQL::DBAdaptor');

  my $added_alias = 0;
  if(! $self->alias_exists($species)) {
    $self->add_alias($species, $species);
    $added_alias = 1;
  }
  
  if($self->DBAdaptor_exists(@args)) {
    $self->remove_alias($species) if $added_alias;
    throw "We already have a DBAdaptor registered under the species '$species' and '$group'";
  }
  
  $self->{db_adaptor}->{$species}->{$group} = $adaptor;
  return;
}

sub add_DNAAdaptor {
  my ($self, @args) = @_;
  
  my ($adaptor, $species, $group, $dna_adaptor, $dna_species, $dna_group) = 
    $self->_process_args([qw/dna_adaptor dna_species dna_group/], @args);
  ($dna_adaptor, $dna_species, $dna_group) = 
    $self->_process_db_adaptor($dna_adaptor, $dna_species, $dna_group);
  
  throw "No species given" if ! defined $species;
  throw "No group given" if ! defined $group;
  throw "No dna species given" if ! defined $dna_species;
  throw "No dna group given" if ! defined $dna_species;
  
  $self->{dna_adaptor}->{$species}->{$group} = [$dna_species, $dna_group];
  
  return;
}

sub add_alias {
  my ($self, $species, $alias) = @_;
  $alias = wrap_array($alias);
  my $lc_species = lc($species);
  foreach my $a (map { lc($_) } @{$alias}) {
    $self->{alias}->{$a} = $lc_species;
  } 
  return;
}

### Removal methods

sub remove_DBAdaptor {
  my ($self, @args) = @_;
  my ($adaptor, $species, $group) = $self->_process_args([], @args);
  throw "No species given" if ! defined $species;
  throw "No group given" if ! defined $group;
  delete $self->{db_adaptor}->{$species}->{$group};
  return;
}

sub remove_alias {
  my ($self, $alias) = @_;
  delete $self->{alias}->{lc($alias)};
  return;
}

### Merging

sub merge {
  my ($self, $other_registry, $verbose) = @_;
  
  assert_ref($other_registry, 'Bio::EnsEMBL::Registry::Storage');
  
  my $no_alias_lookup = 1;
  foreach my $db_adaptor (@{$other_registry->get_all_DBAdaptors()}) {
    my ($species, $group) = $self->_process_db_adaptor($db_adaptor, $no_alias_lookup);
    if(exists $self->{db_adaptor}->{$species}->{$group}) {
      if($verbose) {
        printf("Species '%s' and group '%s' is already stored\n", $species, $group);
      }
    }
    else {
      $self->{db_adaptor}->{$species}->{$group} = $db_adaptor;
    }
  }
  
  foreach my $alias (keys %{$other_registry->{alias}}) {
    if(! $self->alias_exist($alias)) {
      my $species = $other_registry->{alias}->{$alias};
      $self->add_alias($species, $alias);
    }
  }
  
  return;
}

### Generic methods

sub _process_args {
  my ($self, $extra_args, @args) = @_;
  my ($adaptor, $species, $group, @extras) = rearrange([qw/adaptor species group/, @{$extra_args}], @args);
  ($adaptor, $species, $group) = $self->_process_db_adaptor($adaptor, $species, $group);
  return ($adaptor, $species, $group, @extras);
}

sub _process_db_adaptor {
  my ($self, $adaptor, $species, $group, $no_alias_lookup) = @_;
  if(check_ref($adaptor, 'Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    $species  ||= $adaptor->species();
    $group    ||= $adaptor->group();
  }
  if(!$no_alias_lookup) {
    my $decoded_alias = $self->get_alias($species);
    $species = $decoded_alias if defined $decoded_alias;
  }
  return ($adaptor, $species, $group);
}

1; 