package Bio::EnsEMBL::Registry::Loader::Base;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(check_ref);

my $CORE_DBADAPTOR = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $REGULATION_DBADAPTOR = 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor';

sub new {
  my ($class, $registry, $verbose, $no_caching) = @_;
  throw "No registry defined" unless $registry;
  
  $class = ref($class) || $class;
  my $self = bless({ 
    imported_packages => {}, 
    verbose => $verbose, 
    registry => $registry,
    no_caching => $no_caching }, $class);
  
  return $self;
}

sub registry {
  my ($self, $registry) = @_;
  $self->{registry} = $registry if defined $registry;
  return $self->{registry};
}

sub verbose {
  my ($self, $verbose) = @_;
  $self->{verbose} = $verbose if defined $verbose;
  return $self->{verbose};
}

sub no_caching {
  my ($self, $no_caching) = @_;
  $self->{no_caching} = $no_caching if defined $no_caching;
  return $self->{no_caching};
}

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
    userupload    => $CORE_DBADAPTOR,
    variation     => 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
    vega          => $CORE_DBADAPTOR,
  }->{$group};
}

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