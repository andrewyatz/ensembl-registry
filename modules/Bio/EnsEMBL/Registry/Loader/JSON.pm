package Bio::EnsEMBL::Registry::Loader::JSON;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Registry::Loader::PerlHash/;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use JSON;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  #Setup the JSON serialiser object
  $self->{json} = JSON->new();
  $self->{json}->relaxed();
  $self->{json}->pretty();
  $self->{json}->canonical();
  
  #Return outself
  return $self;
}

=head2 load_registry

  Arg [1]    : String $json raw JSON document as a scalar
  Example    : my $json = slurp('/path/to/json'); # bring the JSON document into memory
               $json_loader->load_registry($json);
  Description: Parses the given JSON document and loads the
               associcated Registry with the specified adaptors
  Returntype : None
  Exceptions : Thrown if $json was not defined or it was empty
  Caller     : General
  Status     : At risk

=cut

sub load_registry {
  my ($self, $json) = @_;
  throw "JSON was not defined" if ! defined $json;
  throw "JSON was empty" if ! $json;
  my $perl_hash = $self->{json}->decode($json);
  return $self->SUPER::load_registry($perl_hash);
}

=head2 serialise_registry

  Example    : my $json_doc = $json_loader->serialise_registry();
  Description: For the associated registry object this methods returns the String
               JSON document representation
  Returntype : String of the serialised registry
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub serialise_registry {
  my ($self) = @_;
  my $perl_hash = $self->SUPER::serialise_registry();  
  return $self->{json}->encode($perl_hash);
}

1;