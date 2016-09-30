package CSG::Mapper::Command::import;

use CSG::Mapper -command;
use CSG::Base qw(parsers file);
use CSG::Constants qw(:basic :mapping);
use CSG::Mapper::Config;
use CSG::Mapper::DB;

Readonly::Array my @IMPORT_FIELDS => (
  qw(
    center
    run_dir
    filename
    study
    pi
    sample_id
    state_b37
    state_b38
    year
    flagstat
    fullpath
    )
);

sub opt_spec {
  ## no tidy
  return (
    ['filename|f=s', 'Filename to import or - to read from stdin'],
    ['headers',      'Import file contains a header row'         ],
  );
  ## use tidy
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{filename} ne $DASH and not -e $opts->{filename}) {
    $self->usage_error('Unable to locate import filename on disk');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $config  = CSG::Mapper::Config->new(project => $self->app->global_options->{project});
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
    my $path     = abs_path($line->fullpath);
    my $file     = Path::Class->file($path);
    my @comps    = $file->components();
    my $hostname = $comps[3] // $project;

    my $proj = $schema->resultset('Project')->find_or_create({name => $project});
    my $center = $schema->resultset('Center')->find_or_create({name => $line->center});
    my $study = $schema->resultset('Study')->find_or_create({name => $line->study});
    my $host = $schema->resultset('Host')->find_or_create({name => $hostname});
    my $pi = $schema->resultset('Pi')->find_or_create({name => $line->pi});

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
        year       => $line->year,
        flagstat   => $line->flagstat,
      }
    );
  }
}

1;

__END__

=head1

CSG::Mapper::Command::import - import remapping jobs
