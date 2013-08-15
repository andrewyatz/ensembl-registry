package Bio::EnsEMBL::RegistryAdditions;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::StableIdsDBAadptorLookup;

#################################################################################
###### None of these methods are meant for this class but for the Registry ######
#################################################################################

=head2 get_species_and_object_type

  Description:  Get the species name, object type (gene, transcript,
                translation, or exon etc.), and database type for a
                stable ID.

  Arg [1]    :  String stable_id
                The stable ID to find species and object type for.

  Arg [2]    :  String known_type (optional)
                The object type of the stable ID, if it is known e.g. gene

  Arg [3]    :  String known_species (optional)
                The species, if known e.g. human

  Arg [4]    :  String known_group (optional)
                The database group, if known
                
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
  my $scanner = Bio::EnsEMBL::Utils::StableIdsDBAadaptorLookup->new($self);
  return $scanner->fetch_object_location($stable_id, $known_type, $known_species, $known_group);
}

1;