package CSG::Mapper::Command::stat;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Mapper::DB;

use Text::ASCIITable;

my $schema  = CSG::Mapper::DB->new();

sub opt_spec {
  return (
    ['job-id=s',       'job to provide stats for'],
    ['time-left',      'calculate time remaining in hours for a given jobid'],
    ['totals',         'various counts'],
    ['step|s=s',       'step a result is in'],
    ['avg-job-time=s', 'average job run time for a given step']
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
    $self->_totals($self->app->global_options->{build});
  }

  if ($opts->{avg_job_time}) {
    $self->_avg_job($opts->{avg_job_time});
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

  my $project = $schema->resultset('Project')->find({name => $self->app->global_options->{project}});

  my $table   = Text::ASCIITable->new({headingText => 'Total Samples: ' . $project->samples->count});
  my @columns = (qw(requested submitted started completed cancelled failed)); 
  $table->setCols(('Step', map {ucfirst($_)} @columns));

  for my $step ($schema->resultset('Step')->all) {
    my $results = $schema->resultset('ResultsStatesStep')->current_results_by_step($build, $step->name);

    my %totals = map {$_ => 0} @columns;
    $totals{$_->state->name}++ for $results->all;

    $table->addRow($step->name, @totals{@columns});
  }

  $table->addRowLine();

  print $table;
}

sub _avg_job {
  my ($self, $step) = @_;

  my $time    = 0;
  my $build   = $self->app->global_options->{build};
  my $results = $schema->resultset('ResultsStatesStep')->current_results_by_step_state($build, $step, 'completed');

  for my $result ($results->all) {
    my $started  = $result->job->started_at;
    my $ended    = $result->job->ended_at;

    next unless $started and $ended;
    my $duration = $ended - $started;

    $time += $duration->in_units('hours');
  }

  printf "Average job run time for %s is: %.2f hrs\n", $step, ($time / $results->count);
}

1;

__END__

=head1

CSG::Mapper::Command::stat - stat remapping jobs
