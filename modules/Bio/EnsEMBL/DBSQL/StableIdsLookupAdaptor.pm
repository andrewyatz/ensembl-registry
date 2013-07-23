package Bio::EnsEMBL::DBSQL::StableIdsLookupAdaptor;

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

Bio::EnsEMBL::DBSQL::StableIdsLookupAdaptor

=head1 DESCRIPTION

Provides a single lookup method which can perform searches of the stable_ids database this
adaptor was constructed against.

=cut

use strict;
use warnings;
use base qw/Bio::EnsEMBL::DBSQL::BaseAdaptor/;

use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 fetch_object_location

  Arg[1]      : String stable identifier to lookup in the current stable_ids DB
  Arg[2]      : (optional) String Specify the known object type e.g. gene, transcript
  Arg[3]      : (optional) String Specify the known species. Automatic alias lookup performed if specified e.g. homo_sapiens, human
  Arg[4]      : (optional) String Specify the known database group e.g. core, vega
  Description : Performs a lookup of the stable_ids schema this adaptor was constructed
                against to locate your specified ID into a species, database type and object type.
                This information can be used in conjunction with the registry to materialise
                the Ensembl object instance.
  Example     : my $id = 'ENSG000000001';
                my ($species, $type, $group) = $stla->fetch_object_location($id);
                my $adaptor = $reg->get_adaptor($species, $group, $type);
                my $object = $adaptor->fetch_by_stable_id($id);
                
                #Or if you need to specify more information
                my ($species, $type, $group) = $stla->fetch_object_location('ARANDOMGENEID', 'gene', 'human', 'otherfeatures');
                
  Exceptions  : Thrown in the event of an unspecified stable id or due to DBI errors
  Returntype  : List of the found species, object type and database group

=cut

sub fetch_object_location {
  my ($self, $stable_id, $known_type, $known_species, $known_db_type) = @_;
  throw "No stable id given" unless $stable_id;
  
  my $sql = 'SELECT name, object_type, db_type FROM stable_id_lookup join species using(species_id) WHERE stable_id = ?';
  my @params = ($stable_id);
  if ($known_species) {
    $sql .= ' AND name = ?';
    push(@params, $known_species);
  }
  if ($known_db_type) {
    $sql .= ' AND db_type = ?';
    push(@params, $known_db_type);
  }
  if ($known_type) {
    $sql .= ' AND object_type = ?';
    push(@params, $known_type);
  }
  
  my $results = $self->dbc->sql_helper->execute(-SQL => $sql, -PARAMS => \@params)
  return unless @{$results};
  my ($species, $type, $db_type) = @{$results->[0]};
  return ($species ,$type, $db_type);
}

1;