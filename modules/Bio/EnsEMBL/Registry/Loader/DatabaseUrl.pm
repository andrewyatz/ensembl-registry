package Bio::EnsEMBL::Registry::Loader::DatabaseUrl;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Registry::Loader::Database/;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::URI qw(parse_uri);

=head2 load_registry

  Arg [1] : string $url The URL of the database server to connect to
  Example : $database_url_loader->load_registry('mysql://anonymous@ensembldb.ensembl.org:3306');
            
            $database_url_loader->load_registry('mysql://anonymous@ensembldb.ensembl.org:3306/72');
            
            $database_url_loader->load_registry(
              'mysql://anonymous@ensembldb.ensembl.org:3306/homo_sapiens_core_65_37?group=core&species=homo_sapiens'
            );
            
            $database_url_loader->load_registry(
              'mysql://anonymous@ensembldb.ensembl.org:3306/homo_sapiens_core_65_37?group=core'
            );
            

  Description: Will load the correct versions of the ensembl
               databases for the software release it can find on
               a database instance into the registry. Also adds
               a set of standard aliases. The url format is:
               mysql://[[username][:password]@]hostname[:port].  You
               can also request a specific version for the databases
               by adding a slash and the version number but your
               script may crash as the API version won't match the
               DB version.
               
               You can also specify a database name which will cause the 
               loading of a single DBAdaptor instance. Parameters are
               mapped from a normal URL parameter set to their DBAdaptor
               equivalent. Group must be defined.
  Exceptions : Thrown if the given URL does not parse according to the above 
               scheme and if the specified database cannot be connected to 
               (see L<Bio::EnsEMBL::Registry::Loader::Database::load_registry()> 
               for more information)
  Status     : Stable
 
=cut

sub load_registry {
  my ( $self, $url ) = @_;
  throw "URL was not defined" if ! defined $url;
  throw "URL was empty" if ! $url;
  
  #Detecting if we are loading an entire registry from a single DB instance
  if ( $url =~ /^mysql\:\/\/  #look for the scheme mysql
                ([^\@]+\@)?   #find the username and password which is anything that is not an at symbol. Optional
                ([^\:\/]+)    #look for host name which is anything not a colon
                (\:\d+)?      #find the port. Optional
                (\/\d+)?      #release number to look for. Optional
                \/?$          #Final optional trailing slashes parameter. Regex now finishes
                /xms ) {
    my $user_pass = $1;
    my $host      = $2;
    my $port      = $3;
    my $version   = $4;

    $user_pass =~ s/\@$//;
    my ( $user, $pass ) = $user_pass =~ m/([^\:]+)(\:.+)?/x; #Username is anything not a colon (:). Pass is the remainder 
    $pass    =~ s/^\://x if ($pass); #remove the : 
    $port    =~ s/^\://x if ($port); #remove the :
    $version =~ s/^\///x if ($version); #remove any trailing /
    
    return $self->SUPER::load_registry(
      -HOST       => $host,
      -USER       => $user,
      -PASS       => $pass,
      -PORT       => $port,
      -DB_VERSION => $version
    );
  }
  
  #If we get here then we have been told to load a single adaptor so do so
  my $verbose = $self->verbose();
  my $no_cache = $self->no_cache();
  
  my $uri = parse_uri($url);
  if($uri && $uri->scheme() eq 'mysql') {
    my %params = $uri->generate_dbsql_params();
    if($params{-DBNAME}) {
      $params{-SPECIES} = $params{-DBNAME} unless $params{-SPECIES};
      $params{-NO_CACHE} = 1 if $no_cache;
      my $group = $params{-GROUP};
      
      my $module = $self->group_to_adaptor($group);
      my $success = $self->import_package($module, 0);
      if(! $success) {
        printf("%s module not found on PERL5LIB so %s databases will be ignored\n", $module, $group) if $verbose;
        next;
      }

      if($verbose) {
        printf("Loading database '%s' from group '%s' with DBAdaptor class '%s' from url %s\n", $params{-DBNAME}, $group, $module, $url);
      }
      $module->new(%params, -NO_CACHE => $no_cache);
      return;
    }
    throw "No dbname was specified in the URL '${url}'";
  }
  throw("Only MySQL URLs are accepted. Given URL was '${url}'");
}

1;