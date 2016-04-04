package CSG::Mapper::Command::stat;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Mapper::DB;

sub opt_spec {
  return (
    ['job-id=s',  'job to provide stats for'],
    ['time-left', 'calculate time remaining in hours for a given jobid'],
    ['totals',    'various counts'],
    ['build|b=s', 'reference build'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{time_left} and not $opts->{job_id} and not $self->app->global_options->{cluster}) {
    $self->usage_error('cluster and job-id are required for the time-left stat');
  }

  if ($opts->{totals} and not $opts->{build}) {
    $self->usage_error('build is required when viewing totals');
  }

  unless ($opts->{build} =~ /37|38/) {
    $self->usage_error('invalid reference build');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  if ($opts->{time_left}) {
    $self->_time_left($opts->{job_id});
  }

  if ($opts->{totals}) {
    $self->_totals($opts->{build});
  }
}

sub _time_left {
  my ($self, $job_id) = @_;

  my $job = CSG::Mapper::Job->new(
    cluster => $self->app->global_options->{cluster},
    job_id  => $job_id,
  );

  say $job->time_remaining();
}

sub _totals {
  my ($self, $build) = @_;

  my $schema  = CSG::Mapper::DB->new();
  my $project = $schema->resultset('Project')->find({name => $self->app->global_options->{project}});

  my $total     = $project->samples->count;
  my $completed = $schema->resultset('Result')->completed($build)->count;
  my $failed    = $schema->resultset('Result')->failed($build)->count;
  my $running   = $schema->resultset('Result')->started($build)->count;
  my $submitted = $schema->resultset('Result')->submitted($build)->count;
  my $cancelled = $schema->resultset('Result')->cancelled($build)->count;
  my $requested = $schema->resultset('Result')->requested($build)->count;

  print << "EOF"
Completed:  $completed
Submitted:  $submitted
Requested:  $requested
Cancelled:  $cancelled
Failed:     $failed
----------
Total:      $total
EOF
}

1;

__END__

=head1

CSG::Mapper::Command::stat - stat remapping jobs
