package Bio::EnsEMBL::Utils::StableIdsDBAdaptorLookup;

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

Bio::EnsEMBL::Utils::StableIdsDBAdaptorLookup

=head1 DESCRIPTION

Provides a mechanism to perform linear scans and queries of databases
in order to locate a stable identifier. This class should not be
used directly. Instead use the lookup methods in 
Bio::EnsEMBL::Registry.

=cut

use strict;
use warnings;

use Scalar::Util qw/looks_like_number/;
use Bio::EnsEMBL::Utils::Exception qw(throw);

=head2 new

  Arg[1]      : Bio::EnsEMBL::Registry The registry to lookup against
  Description : Constructs an instance of the lookup

=cut

sub new {
  my ($class, $registry) = @_;
  throw "Was not given a Registry during construction" if ! $registry;
  $class = ref($class) || $class;
  return bless({ registry => $registry }, $class);
}

=head2 fetch_object_location

  Arg[1]      : String stable identifier to find in the currently available set of DBAdaptors
  Arg[2]      : (optional) String Specify the known object type e.g. gene, transcript
  Arg[3]      : (optional) String Specify the known species. Automatic alias lookup performed if specified e.g. homo_sapiens, human
  Arg[4]      : (optional) String Specify the known database group e.g. core, vega. Defaults to core
  Description : Performs a linear scan of all available DBAdaptors to find the location of
                stable ID in an Ensembl database. Currently we support
                Genes, Transcripts, Exons, Translations, Operons and OperonTranscripts
                from core-like databases. Compara GeneTree identifiers are also supported.
                
                By default we only look in core databases so compara DBs must be forced.
  Example     : my $id = 'ENSG000000001';
                my ($species, $type, $group) = $stla->fetch_object_location($id, undef, undef, 'core');
                my $adaptor = $registry->get_adaptor($species, $group, $type);
                my $object = $adaptor->fetch_by_stable_id($id);
                
                #Or if you need to specify more information
                my ($species, $type, $group) = $stla->fetch_object_location('ARANDOMGENEID', 'gene', 'human', 'otherfeatures');
                
  Exceptions  : Thrown in the event of an unspecified stable id or due to DBI errors
  Returntype  : List of the found species, object type and database group

=cut

sub fetch_object_location {
  my ($self, $stable_id, $known_type, $known_species, $known_db_type) = @_;
  throw "No stable id given" unless $stable_id;
  
  $known_db_type ||= 'core';
  
  my ($types, $statements_lookup);
  if($known_db_type eq 'compara') {
    #If it was a compara group then we only have one type to work with
    $types = defined $known_type ? [$known_type] : ['genetree'];
    $statements_lookup = $self->_compara_statements();
  }
  else {
    #otherwise it is a core DB so we can use many types
    $types = defined $known_type ? [$known_type] : ['gene', 'transcript', 'translation', 'exon', 'operon', 'operontranscript'];
    $statements_lookup = $self->_core_statements();
  }
  
  #Get adaptors
  my %get_dbadaptor_args = ('-GROUP' => $known_db_type);
  $get_dbadaptor_args{'-SPECIES'} = $known_species if $known_species;
  my @dbas = @{$self->{registry}->get_all_DBAdaptors(%get_adaptors_args)};
  
  #Processed DBs are a multi-species only speed-up.
  #Once we've queried a DB we do not need to query it again so store the DBC locator string
  my %processed_db;
  foreach my $dba (@dbas) {
    my $is_multispecies = $dba->is_multispecies();
    if($is_multispecies) {
      my $locator = $dba->dbc->locator;
      next if exists $processed_db{$locator};
    }
    my ($species, $type, $group) = $self->_db_scan($stable_id, $types, $dba, $statements_lookup);
    if($is_multispecies) {
      $processed_db{$dba->dbc->locator()} = 1;
    }
    return ($species, $type, $group) if $species;
  }
  
  return;
}

=head2 _db_scan

  Arg[1]      : String stable id to search with
  Arg[2]      : ArrayRef of known types of objects to lookup for
  Arg[3]      : Bio::EnsEMBL::DBSQL::DBAdaptor to search against
  Arg[4]      : HashRef of type to SQL statement. Must have one sprintf substitution 
                available which is the database name
  Description : Provides a generic database lookup mechanism for stable id finding
  Returntype  : List of the species, object type and adaptor group it was first found in

=cut

sub _db_scan {
  my ($self, $stable_id, $known_types, $dba, $statement_lookup) = @_;
  my $helper = $self->db->dbc->sql_helper();
  foreach my $type (@t{$known_types}) {
    #Find the statement & sprintf it
    my $statement = sprintf($statement_lookup{lc $type}, $dba->dbc->dbname);
    my $species = $helper->execute_simple(
      -SQL => $sql, -PARAMS => [$stable_id], -NO_ERROR => 1
    );
    if($species) {
      #If it was a number then it must have been a simple boolean return 
      #from the statement so make sure we replace it with the right species name
      $species = $dba->species() if looks_like_number($species);
      return ($species, $type, $dba->group());
    }
  }
  return;
}

sub _core_statements {
  return {
    gene => 'SELECT m.meta_value '
      . 'FROM %1$s.gene '
      . 'JOIN %1$s.seq_region USING (seq_region_id) '
      . 'JOIN %1$s.coord_system USING (coord_system_id) '
      . 'JOIN %1$s.meta m USING (species_id) '
      . 'WHERE stable_id = ? '
      . 'AND m.meta_key = "species.production_name"',
    transcript => 'SELECT m.meta_value '
      . 'FROM %1$s.transcript '
      . 'JOIN %1$s.seq_region USING (seq_region_id) '
      . 'JOIN %1$s.coord_system USING (coord_system_id) '
      . 'JOIN %1$s.meta m USING (species_id) '
      . 'WHERE stable_id = ? '
      . 'AND m.meta_key = "species.production_name"',
    exon => 'SELECT m.meta_value '
      . 'FROM %1$s.exon '
      . 'JOIN %1$s.seq_region USING (seq_region_id) '
      . 'JOIN %1$s.coord_system USING (coord_system_id) '
      . 'JOIN %1$s.meta m USING (species_id) '
      . 'WHERE stable_id = ? '
      . 'AND m.meta_key = "species.production_name"',
    translation => 'SELECT m.meta_value '
      . 'FROM %1$s.translation tl '
      . 'JOIN %1$s.transcript USING (transcript_id) '
      . 'JOIN %1$s.seq_region USING (seq_region_id) '
      . 'JOIN %1$s.coord_system USING (coord_system_id) '
      . 'JOIN %1$s.meta m USING (species_id) '
      . 'WHERE tl.stable_id = ? '
      . 'AND m.meta_key = "species.production_name"',
    operon => 'SELECT m.meta_value '
      . 'FROM %1$s.operon '
      . 'JOIN %1$s.seq_region USING (seq_region_id) '
      . 'JOIN %1$s.coord_system USING (coord_system_id) '
      . 'JOIN %1$s.meta m USING (species_id) '
      . 'WHERE stable_id = ? '
      . 'AND m.meta_key = "species.production_name"',
    operontranscript => 'SELECT m.meta_value '
      . 'FROM %1$s.operon_transcript '
      . 'JOIN %1$s.seq_region USING (seq_region_id) '
      . 'JOIN %1$s.coord_system USING (coord_system_id) '
      . 'JOIN %1$s.meta m USING (species_id) '
      . 'WHERE stable_id = ? '
      . 'AND m.meta_key = "species.production_name"',
  };
}

sub _compara_statements {
  return {
    genetree => 'SELECT 1 FROM %1$s.gene_tree_root WHERE stable_id =?',
  };
}

1;