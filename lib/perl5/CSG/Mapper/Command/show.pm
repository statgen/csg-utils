package CSG::Mapper::Command::show;

use CSG::Mapper -command;
use CSG::Base qw(formats);
use CSG::Constants;
use CSG::Mapper::DB;
use CSG::Mapper::Job;
use CSG::Mapper::Sample;

my $schema = CSG::Mapper::DB->new();

sub opt_spec {
  return (
    ['info',          'combined output of job and sample(for backward compatiblity)'],    # TODO - remove after old jobs clear queue
    ['meta-id=i',     'job id for the meta data record (for backward compatibility)'],    # TODO - remove after old jobs clear queue
    ['job-info=i',    'display basic job info'],
    ['sample-info=s', 'display all info on a given sample'],
    ['result-info=s', 'display sample and result info for a given sample id'],
    ['step=s',        'display results for a given step (e.g. bam2fastq, align)'],
    ['state=s', 'display results for a given state (e.g. submitted, requested, failed)'],
    ['stale',   'find any jobs that are no longer queued but still in a running state (i.e. started, submitted)'],
    [
      'format=s',
      'output format (valid format: yaml|txt) [default: yaml]', {
        default   => 'yaml',
        callbacks => {
          regex => sub {
            shift =~ /yaml|txt/;
          }
        }
      }
    ]
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{state}) {
    my $state = $schema->resultset('State')->find({name => $opts->{state}});
    unless ($state) {
      $self->usage_error('invalid state');
    }

    $self->{stash}->{state} = $state;
  }

  if ($opts->{step}) {
    my $step = $schema->resultset('Step')->find({name => $opts->{step}});
    unless ($step) {
      $self->usage_error('invalid step');
    }

    $self->{stash}->{step} = $step;
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  if ($opts->{info}) {
    my $job = $schema->resultset('Job')->find($opts->{meta_id});
    $self->_dump($opts->{format}, $self->_job_info($job));

    my $sample = $job->result->sample;
    $self->_dump($opts->{format}, $self->_sample_info($sample));
  }

  if ($opts->{job_info}) {
    my $job  = $schema->resultset('Job')->find($opts->{job_info});
    my $info = $self->_job_info($job);

    $self->_dump($opts->{format}, $info);
  }

  if ($opts->{sample_info}) {
    my $sample = $schema->resultset('Sample')->find({sample_id => $opts->{sample_info}});
    my $info = $self->_sample_info($sample);

    $self->_dump($opts->{format}, $info);
  }

  if ($opts->{result_info}) {
    my $sample = $schema->resultset('Sample')->find({sample_id => $opts->{result_info}});

    my $info = $self->_sample_info($sample);
    $info->{sample}->{results} = $self->_result_info($sample);

    $self->_dump($opts->{format}, $info);
  }

  if ($opts->{stale}) {
    return $self->_stale();
  }

  if ($opts->{state}) {
    my $build = $self->app->global_options->{build};
    my $state = $self->{stash}->{state};
    my $step  = $self->{stash}->{step};

    for my $result ($schema->resultset('ResultsStatesStep')->current_results_by_step_state($build, $step->name, $state->name)) {
      say $result->result->status_line();
    }
  }
}

sub _dump {
  my ($self, $format, $data) = @_;

  my $formats = {
    txt  => sub {print Dumper shift},
    yaml => sub {print Dump shift},
  };

  croak 'invalid format' unless exists $formats->{$format};

  $formats->{$format}->($data);

  return;
}

sub _sample_info {
  my ($self, $sample) = @_;

  my $result     = $sample->result_for_build($self->app->global_options->{build});
  my $sample_obj = CSG::Mapper::Sample->new(
    cluster => $self->app->global_options->{cluster},
    record  => $sample,
    build   => $result->build,
  );

  return {
    sample => {
      id            => $sample->id,
      sample_id     => $sample->sample_id,
      center        => $sample->center->name,
      study         => $sample->study->name,
      pi            => $sample->pi->name,
      host          => $sample->host->name,
      filename      => $sample->filename,
      run_dir       => $sample->run_dir,
      fullpath      => $sample->fullpath,
      out_dir       => $sample_obj->result_path,
      current_state => $result->current_state,
      current_step  => $result->current_step,
      build         => $result->build,
      fastqs        => [map +{read_group => $_->read_group, path => $_->path}, $sample->fastqs],
    }
  };
}

sub _result_info {
  my ($self, $sample) = @_;

  my @results = ();
  my $result  = $sample->result_for_build($self->app->global_options->{build});

  for ($result->results_states_steps->all) {
    push @results, {
      job_id => $_->job_id,
      state  => $_->state->name,
      step   => $_->step->name,
      };
  }

  return \@results;
}

sub _job_info {
  my ($self, $job) = @_;

  return {
    job => {
      id        => $job->id,
      job_id    => $job->job_id,
      result_id => $job->result->id,
      cluster   => $job->cluster,
      procs     => $job->procs,
      memory    => $job->memory,
      walltime  => $job->walltime,
      node      => $job->node,
      delay     => $job->delay,
      tmp_dir   => $job->tmp_dir,
      submitted => ($job->submitted_at) ? $job->submitted_at->ymd . $SPACE . $job->submitted_at->hms : $EMPTY,
      started   => ($job->started_at) ? $job->started_at->ymd . $SPACE . $job->started_at->hms : $EMPTY,
      created   => $job->created_at->ymd . $SPACE . $job->created_at->hms,
    }
  };
}

sub _stale {
  my ($self) = @_;

  my $step    = $self->{stash}->{step};
  my $cluster = $self->app->global_options->{cluster};
  my $build   = $self->app->global_options->{build};
  my $results = $schema->resultset('ResultsStatesStep')->current_results_by_step($build, $step->name);

  for my $result ($results->all) {
    next unless $result->state->name eq 'started';
    next unless $result->job->cluster eq $cluster;

    my $job = CSG::Mapper::Job->new(
      cluster => $cluster,
      job_id  => $result->job->job_id
    );

    my $job_state = $job->state;
    next if $job_state eq 'running';
    say $result->result->status_line . 'JOBID: ' . $job->job_id . ' JOBSTATUS: ' . $job_state;
  }
}

1;

__END__

=head1

CSG::Mapper::Command::show - show remapping jobs
