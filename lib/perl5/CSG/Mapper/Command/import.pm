package CSG::Mapper::Command::import;

use CSG::Mapper -command;
use CSG::Base qw(parsers file);
use CSG::Constants qw(:basic :mapping);
use CSG::Mapper::Config;
use CSG::Mapper::DB;

Readonly::Array my @IMPORT_FIELDS => (qw(center run_dir filename study pi sample_id state_b37 state_b38 fullpath));

sub opt_spec {
  return (
    ['filename|f=s', 'Filename to import or - to read from stdin'],
    ['headers', 'Import file contains a header row'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($self->app->global_options->{cluster}) {
    $self->usage_error('Cluster environment is required');
  }

  unless ($self->app->global_options->{cluster} =~ /$VALID_CLUSTER_REGEXPS/) {
    $self->usage_error('Invalid cluster environment');
  }

  unless ($self->app->global_options->{project}) {
    $self->usage_error('Project is required');
  }

  my $config = CSG::Mapper::Config->new(project => $self->app->global_options->{project});
  $self->{stash}->{config} = $config;

  if ($opts->{filename} ne $DASH and not -e $opts->{filename}) {
    $self->usage_error('Unable to locate import filename on disk');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $config  = $self->{stash}->{config};
  my $schema  = CSG::Mapper::DB->new();
  my $project = $self->app->global_options->{project};
  my $cluster = $self->app->global_options->{cluster};

  my %parse_params = (
    filename => $opts->{filename},
    fields   => \@IMPORT_FIELDS,
  );

  if ($opts->{filename} eq $DASH) {
    delete $parse_params{filename};
    $parse_params{filehandle} = io->stdin->tie;
  }

  my $csv;
  try {
    $csv = Class::CSV->parse(%parse_params);
  }
  catch {
    croak "failed to parse csv: $_";
  };

  my @lines = @{$csv->lines()};
  shift @lines if $opts->{headers};

  for my $line (@lines) {
    my $hostname    = $project;
    my $incoming    = File::Spec->join($config->get($cluster, 'prefix'), $project, $config->get($project, 'incoming_dir'));
    my $center_path = File::Spec->join($incoming, $line->center);

    if (-l $center_path) {
      my $file  = Path::Class->file(readlink($center_path));
      my @comps = $file->components();

      $hostname = $comps[4];
    }

    my $proj   = $schema->resultset('Project')->find_or_create({name => $project});
    my $center = $schema->resultset('Center')->find_or_create({name => $line->center});
    my $study  = $schema->resultset('Study')->find_or_create({name => $line->study});
    my $host   = $schema->resultset('Host')->find_or_create({name => $hostname});
    my $pi     = $schema->resultset('Pi')->find_or_create({name => $line->pi});

    $schema->resultset('Sample')->find_or_create(
      {
        sample_id  => $line->sample_id,
        center_id  => $center->id,
        study_id   => $study->id,
        pi_id      => $pi->id,
        host_id    => $host->id,
        project_id => $proj->id,
        filename   => $line->filename,
        run_dir    => $line->run_dir,
        fullpath   => $line->fullpath,
      }
    );
  }
}

1;

__END__

=head1

CSG::Mapper::Command::import - import remapping jobs
