package Bio::EnsEMBL::Registry::Loader::ConfigIni;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Registry::Loader::PerlHash/;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

#Do IniFile parser detection. Lack of it means we cannot load
#the config file
our $INIFILE_OK = 0;
eval { 
  require Config::IniFiles; 
  $INIFILE_OK = 1;
};

=head2 load_registry

  Arg[1]     : String path to the ini file. Always assumed to be local. Behaves
               according to the rules in Config::IniFiles (so it can be a path,
               *GLOB reference, *GLOB or IO::File object)
  Example    : $config_ini->load_registry('/path/to/config.ini');
  Description: Parses the given path into a Config::IniFiles object and loads
               all adaptors into the associcated Registry
  Returntype : None
  Exceptions : Thrown if Config::IniFiles is not available on PERL5LIB
  Caller     : General
  Status     : At risk

=cut

sub load_registry {
  my ($self, $path) = @_;
  throw "Cannot load from an IniFile since Config::IniFiles is not available on PERL5LIB" if ! $INIFILE_OK;
  my $cfg = Config::IniFiles->new(-file => $path);
  return unless $cfg;
  return $self->_load_from_config_ini($cfg);
}

sub _load_from_config_ini {
  my ($self, $cfg) = @_;
  
  my $registry= $self->registry();
  my $verbose = $self->verbose();
  my $no_caching = $self->no_caching();
  
  my %aliases;
  my @dba_hashes;
  my %default_adaptor_args;
  
  #TODO Config::IniFiles has a "-default" construction option. Investigate
  if ($cfg->SectionExists('default')) {
    my $alias = $cfg->val( 'default', 'alias' );
    $cfg->delval( 'default', 'alias' );
    my $species = $cfg->val( 'default', 'species' );
    if ( defined $alias && defined $species ) {
      $aliases{$species} = [split(/\r?\n/, $alias)];
    }
    %default_adaptor_args = map { $_, $cfg->val('default', $_) } $cfg->Parameters('default'); 
  }
  
  foreach my $section ( $cfg->Sections() ) {
    next if $section eq 'default';
    my $group = $cfg->val($section, 'group') || $default_adaptor_args{'group'};
    if (!defined($group)) {
      warning sprintf("Key 'group' is undefined for configuration section '%s', skipping this section.\n", $section);
      next;
    }
    
    my $species = $cfg->val($section, 'species');
    if(!defined($species)) {
      warning sprintf("Key 'species' is undefined for configuration section '%s', skipping this section.\n", $section);
      next;
    }
    
    my $alias = $cfg->val($section, 'alias');
    $cfg->delval($section, 'alias');
    
    my %adaptor_args = %default_adaptor_args;
    foreach my $parameter ( $cfg->Parameters($section) ) {
      $adaptor_args{$parameter} = $cfg->val($section, $parameter);
    }
    my $split_aliases = ($alias) ? [split(/\r?\n/, $alias)] : [];
    if(! exists $aliases{$species}) {
      $aliases{$species} = [];
    }
    push(@{$aliases{$species}}, @{$split_aliases});

    push(@dba_hashes, {%adaptor_args, -NO_CACHE => $no_caching});
  }
  
  $self->SUPER::load_registry({ adaptors => \@dba_hashes, aliases => \%aliases});
  
  return;
}

=head2 serialise_registry

  Example    : my $config_ini_obj = $config_ini->serialise_registry();
  Description: For the associated registry object this methods returns the ConfigIni
               representation ready for serialisation to disk
  Returntype : Config::IniFiles object
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub serialise_registry {
  my ($self) = @_;
  my $registry = $self->registry();
  my $cfg = Config::IniFiles->new();
  my @dbas = sort { $a->species() cmp $b->species() } @{$registry->get_all_DBAdaptors()};
  foreach my $dba (@dbas) {
    my $section = sprintf('%s_%s', $dba->species(), $dba->group());
    $cfg->AddSection($section);
    $cfg->newval($section, 'species', $dba->species());
    $cfg->newval($section, 'species_id', $dba->species_id());
    $cfg->newval($section, 'group', $dba->group());
    $cfg->newval($section, 'multispecies_db', ($dba->is_multispecies() ? 1 : 0) );
    
    my $dbc = $dba->dbc();
    $cfg->newval($section, 'dbname', $dbc->dbname());
    $cfg->newval($section, 'host', $dbc->host());
    $cfg->newval($section, 'user', $dbc->username());
    $cfg->newval($section, 'pass', $dbc->password()) if $dbc->password();
    $cfg->newval($section, 'port', $dbc->port());
    $cfg->newval($section, 'driver', $dbc->driver());
    $cfg->newval($section, 'disconnect_when_inactive', ($dbc->disconnect_when_inactive() ? 1 : 0));
    $cfg->newval($section, 'wait_timeout', $dbc->timeout());
    $cfg->newval($section, 'reconnect_when_connection_lost', ($dbc->reconnect_when_lost() ? 1 : 0));
    
    my $aliases = $registry->get_all_aliases($dba->species());
    $cfg->newval($section, 'alias', $aliases);
  }
  return $cfg;
}

1;