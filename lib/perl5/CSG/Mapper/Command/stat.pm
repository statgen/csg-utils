package CSG::Mapper::Command::stat;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Mapper::DB;

sub opt_spec {
  return (
    ['job-id=s',  'job to provide stats for'],
    ['time-left', 'calculate time remaining in hours for a given jobid'],
    ['totals',    'various counts'],
    ['step|s=s', 'step a result is in', {required => 1}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{time_left} and not $opts->{job_id}) {
    $self->usage_error('job-id id required for the time-left stat');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  if ($opts->{time_left}) {
    $self->_time_left($opts->{job_id});
  }

  if ($opts->{totals}) {
    $self->_totals($opts->{step}, $self->app->global_options->{build});
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
  my ($self, $step, $build) = @_;

  my $schema  = CSG::Mapper::DB->new();
  my $project = $schema->resultset('Project')->find({name => $self->app->global_options->{project}});

  my $total     = $project->samples->count;
  my $completed = $schema->resultset('Result')->completed($step, $build)->count;
  my $failed    = $schema->resultset('Result')->failed($step, $build)->count;
  my $running   = $schema->resultset('Result')->started($step, $build)->count;
  my $submitted = $schema->resultset('Result')->submitted($step, $build)->count;
  my $cancelled = $schema->resultset('Result')->cancelled($step, $build)->count;
  my $requested = $schema->resultset('Result')->requested($step, $build)->count;

  print << "EOF"
STEP: $step
----------
Requested:  $requested
Submitted:  $submitted
Running:    $running
Completed:  $completed
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
