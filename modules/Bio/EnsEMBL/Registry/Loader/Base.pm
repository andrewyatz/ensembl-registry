package Bio::EnsEMBL::Registry::Loader::Base;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(check_ref);

my $CORE_DBADAPTOR = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $REGULATION_DBADAPTOR = 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor';

=head2 new

  Arg[1]      : Bio::EnsEMBL::Registry The Registry instance to interact with
  Arg[2]      : Boolean Switch verbose statements on or off
  Arg[3]      : Boolean Switch caching on or off
  Description : Builds a Base instance of a Loader object
  Exceptions  : Thrown if registry was not available
  Returntype  : Bio::EnsEMBL::Registry::Loader::Base instance

=cut

sub new {
  my ($class, $registry, $verbose, $no_cache) = @_;
  throw "No registry defined" unless $registry;
  
  $class = ref($class) || $class;
  my $self = bless({ 
    imported_packages => {}, 
    verbose => $verbose, 
    registry => $registry,
    no_cache => $no_cache }, $class);
  
  return $self;
}

=head2 registry

  Arg[1]      : Bio::EnsEMBL::Registry The Registry instance to store
  Description : Getter/setter for the backing registry object
  Returntype  : Bio::EnsEMBL::Registry

=cut

sub registry {
  my ($self, $registry) = @_;
  $self->{registry} = $registry if defined $registry;
  return $self->{registry};
}

=head2 verbose

  Arg[1]      : Boolean The verbose boolean to store
  Description : Getter/setter for the verbose flag
  Returntype  : Boolean

=cut

sub verbose {
  my ($self, $verbose) = @_;
  $self->{verbose} = $verbose if defined $verbose;
  return $self->{verbose};
}

=head2 no_cache

  Arg[1]      : Boolean The no cache boolean to store
  Description : Getter/setter for the no_cache flag. Results in
                the addition of -NO_CACHE to DBAdaptor construction
                calls
  Returntype  : Boolean

=cut

sub no_cache {
  my ($self, $no_cache) = @_;
  $self->{no_cache} = $no_cache if defined $no_cache;
  return $self->{no_cache};
}

=head2 group_to_adaptor

  Arg[1]      : String The group to process
  Description : Holds a look up of group name to DBAdaptor instance which
                should be used to construct if this group is encountered
  Returntype  : String Class name to use

=cut

sub group_to_adaptor {
  my ($self, $group) = @_;
  return {
    ancestral     => $CORE_DBADAPTOR,
    blast         => 'Bio::EnsEMBL::External::BlastAdaptor',
    compara       => 'Bio::EnsEMBL::Compara::DBSQL::DBAdaptor',
    core          => $CORE_DBADAPTOR,
    cdna          => $CORE_DBADAPTOR,
    estgene       => $CORE_DBADAPTOR,
    funcgen       => $REGULATION_DBADAPTOR,
    regulation    => $REGULATION_DBADAPTOR,
    haplotype     => 'Bio::EnsEMBL::ExternalData::Haplotype::DBAdaptor',
    hive          => 'Bio::EnsEMBL::Hive::DBSQL::DBAdaptor',
    ontology      => 'Bio::EnsEMBL::DBSQL::OntologyDBAdaptor',
    otherfeatures => $CORE_DBADAPTOR,
    pipeline      => 'Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor',
    rnaseq        => $CORE_DBADAPTOR,
    snp           => 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor',
    stable_ids    => $CORE_DBADAPTOR,
    userupload    => $CORE_DBADAPTOR,
    variation     => 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
    vega          => $CORE_DBADAPTOR,
  }->{$group};
}

=head2 import_package

  Arg[1]      : String Package to import
  Arg[2]      : Boolean Raise an error if we cannot import the package
  Arg[3]      : Boolean Warn about the inability to import
  Description : Performs a String eval of the given package to import it
                into our %INC space. Results of repeated calls are
                remembered meaning we do not always String eval
  Exceptions  : Thrown if raise error was on
  Returntype  : Boolean indicating success or failure

=cut

sub import_package {
  my ($self, $package, $raise_error, $verbose) = @_;
  throw "Cannot import a package which is not set" if ! $package;
  #Return early if we have already attempted to load the package before
  return $self->{imported_packages}->{$package} if exists $self->{imported_packages}->{$package}; 
  my $eval_result = eval "require $package";
  my $result = 1;
  if(!$eval_result || $@) {
    my $msg = "Error occured whilst requiring the package '$package': $@";
    throw $msg if $raise_error;
    warning $msg if $verbose;
    $result = 0;
  }
  #Remember we have already imported this once
  $self->{imported_packages}->{$package} = $result;
  return $result;
}

1;