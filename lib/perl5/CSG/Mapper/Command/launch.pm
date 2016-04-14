## no critic (NamingConventions::Capitalization, Subroutines::RequireFinalReturn)
package CSG::Mapper::Command::launch;

use CSG::Mapper -command;
use CSG::Base qw(file templates);
use CSG::Constants qw(:basic :mapping);
use CSG::Mapper::Config;
use CSG::Mapper::DB;
use CSG::Mapper::Job;
use CSG::Mapper::Logger;
use CSG::Mapper::Sample;

sub opt_spec {
  return (
    ['limit|l=i',    'Limit number of jobs to submit'],
    ['procs|p=i',    'Number of cores to request'],
    ['memory|m=i',   'Amount of memory to request, in MB'],
    ['walltime|w=i', 'Amount of wallclock time for this job'],
    ['delay=i',      'Amount of time to delay exection in seconds'],
    ['build|b=i',    'Reference build to use (ie; 37 or 38)'],
    ['step=s',       'Job step to launch [bam2fastq|align|all]'],
    ['meta-id=i',    'Job meta record for parent job'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Mapper::DB->new();

  $self->{stash}->{schema} = $schema;

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

  unless ($opts->{step}) {
    $self->usage_error('Step is required');
  }


  my $step = $schema->resultset('Step')->find({name => $opts->{step}});
  $self->{stash}->{step} = $step;
  unless ($step) {
    $self->usage_error('Invalid step specified');
  }

  if ($opts->{meta_id}) {
    my $meta = $schema->resultset('Job')->find($opts->{meta_id});
    unless ($meta) {
      $self->usage_error('Invalid meta job id');
    }
    $self->{stash}->{meta} = $meta;
  }

  if ($opts->{'tmp-dir'}) {
    unless (-e $opts->{'tmp-dir'} and -r $opts->{'tmp-dir'}) {
      $self->usage_error('Temporary disk space does not exist or is not writable');
    }
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $debug   = $self->app->global_options->{debug};
  my $verbose = $self->app->global_options->{verbose};
  my $cluster = $self->app->global_options->{cluster};
  my $project = $self->app->global_options->{project};

  my $jobs   = 0;
  my $schema = $self->{stash}->{schema};
  my $config = $self->{stash}->{config};
  my $step   = $self->{stash}->{step};

  my $project_dir = qq{$FindBin::Bin/..};
  my $prefix      = $config->get($cluster, 'prefix');
  my $workdir     = $config->get($project, 'workdir');

  ## no tidy
  my $procs = ($opts->{procs})
    ? $opts->{procs}
    : ($config->get($cluster, $step->name . '_procs'))
      ? $config->get($cluster, $step->name . '_procs')
      : $config->get($project, 'procs');

  my $memory = ($opts->{memory})
    ? $opts->{memory}
    : ($config->get($cluster, $step->name . '_memory'))
      ? $config->get($cluster, $step->name . '_memory')
      : $config->get($project, 'memory');

  my $walltime = ($opts->{walltime})
    ? $opts->{walltime}
    : ($config->get($cluster, $step->name . '_walltime'))
      ? $config->get($cluster, $step->name . '_walltime')
      : $config->get($project, 'walltime');
  ## use tidy

  my $build   = $opts->{build}    // $config->get($project, 'build');
  my @samples = ();
  my $dep_job_meta = $self->{stash}->{meta};
  if ($dep_job_meta) {
    push @samples, $dep_job_meta->result->sample;
  } else {
    @samples = $schema->resultset('Sample')->search({},{order_by => 'RAND()'});
  }

  for my $sample (@samples) {
    last if $opts->{limit} and $jobs >= $opts->{limit};

    my $logger     = CSG::Mapper::Logger->new();
    my $sample_obj = CSG::Mapper::Sample->new(
      cluster => $cluster,
      record  => $sample,
      build   => $build
    );

    try {
      $sample_obj->incoming_path;
    }
    catch {
      if (not ref $_) {
        $logger->critical('Uncaught exception');
        $logger->debug($_) if $debug;

      } elsif ($_->isa('CSG::Mapper::Exceptions::Sample::NotFound')) {
        $logger->error($_->description);
        $logger->debug('bam_path: ' . $_->bam_path) if $debug;
        $logger->debug('cram_path: ' . $_->cram_path) if $debug;
      } elsif ($_->isa('CSG::Mapper::Exceptions::Sample::SlotFailed')) {
        $logger->error($_->error);
      } else {
        if ($_->isa('Exception::Class')) {
          chomp(my $error = $_->error);
          $logger->critical($error);
        } else {
          $logger->critical('something went sideways');
          print STDERR Dumper $_ if $debug;
        }
      }
    }
    finally {
      unless (@_) {
        $logger->debug('incoming_path: ' . $sample_obj->incoming_path) if $debug;
      }
    };

    next unless $sample_obj->has_incoming_path;

    my $result  = $sample->results->search({build => $build})->first;
    my $tmp_dir = File::Spec->join($config->get($cluster, 'tmp_dir'), $project, $sample_obj->build_str, $sample->sample_id);

    if ($opts->{step} eq 'all') {
      $tmp_dir =  File::Spec->join($config->get($cluster, 'tmp_dir'), $project, $sample_obj->build_str, $opts->{step});
    }

    unless ($dep_job_meta) {
      unless ($result) {
        $result = $sample->add_to_results(
          {
            build    => $build,
            state_id => $schema->resultset('State')->find({name => 'requested'})->id,
          }
        );
      }

      next if $result->state->name ne 'requested';
      next if $result->build ne $build;
    }

    my $delay    = $opts->{delay} // int(rand($MAX_DELAY));
    my $job_meta = $result->add_to_jobs(
      {
        cluster  => $cluster,
        procs    => $procs,
        memory   => $memory,
        walltime => $walltime,
        delay    => $delay,
        step_id  => $step->id,
      }
    );

    $logger->job_id($job_meta->id);

    unless (-e $sample_obj->result_path) {
      $logger->debug('creating out_dir');
      make_path($sample_obj->result_path);
    }

    if ($debug) {
      $logger->debug("cluster: $cluster");
      $logger->debug("procs: $procs");
      $logger->debug("memory: $memory");
      $logger->debug("walltime: $walltime");
      $logger->debug("delay: $delay");
      $logger->debug("build: $build");
      $logger->debug("tmp_dir: $tmp_dir");
      $logger->debug('out_dir: ' . $sample_obj->result_path);
    }

    my $basedir = File::Spec->join($prefix, $workdir);
    $logger->debug("basedir: $basedir") if $debug;
    unless (-e $basedir) {
      make_path($basedir);
      $logger->debug('created basedir') if $debug;
    }

    my $log_dir = $sample_obj->log_dir;
    $logger->debug("log_dir: $log_dir") if $debug;
    unless (-e $log_dir) {
      make_path($log_dir);
      $logger->debug('created log_dir') if $debug;
    }

    my $run_dir = $sample_obj->state_dir;
    $logger->debug("run_dir: $run_dir") if $debug;
    unless (-e $run_dir) {
      make_path($run_dir);
      $logger->debug('created run_dir') if $debug;
    }

    my $gotcloud_conf = File::Spec->join($project_dir, $config->get($cluster, 'gotcloud_conf') . $PERIOD . $sample_obj->build_str);
    $logger->debug("gotcloud conf: $gotcloud_conf") if $debug;
    unless (-e $gotcloud_conf) {
      croak qq{Unable to locate GOTCLOUD_CONF [$gotcloud_conf]};
    }

    my $gotcloud_root = File::Spec->join($basedir, $config->get($cluster, 'gotcloud_root'));
    $logger->debug("gotcloud root: $gotcloud_root") if $debug;
    unless (-e $gotcloud_root) {
      croak qq{GOTCLOUD_ROOT [$gotcloud_root] does not exist!};
    }

    my $gotcloud_ref = File::Spec->join($prefix, $config->get('gotcloud', qq{build${build}_ref_dir}));
    $logger->debug("gotcloud ref_dir: $gotcloud_ref") if $debug;
    unless (-e $gotcloud_ref) {
      croak qq{GOTCLOUD_REF_DIR [$gotcloud_ref] does not exist!};
    }

    my $job_file = File::Spec->join($run_dir, join($DASH, ($sample_obj->sample_id, $step->name, $sample_obj->build_str, $cluster . '.sh')));
    my $tt = Template->new(INCLUDE_PATH => qq($project_dir/templates/batch/$project));

    $tt->process(
      $step->name . q{.sh.tt2}, {
        job => {
          procs      => $procs,
          memory     => $memory,
          walltime   => $walltime,
          build      => $build,
          email      => $config->get($project, 'email'),
          job_name   => join($DASH, ($project, $step->name, $sample_obj->build_str, $sample_obj->sample_id)),
          account    => $config->get($cluster, 'account'),
          workdir    => $log_dir,
          job_dep_id => ($dep_job_meta) ? $dep_job_meta->job_id : undef,
          nodelist   => ($dep_job_meta) ? $dep_job_meta->node   : undef,
        },
        settings => {
          tmp_dir         => $tmp_dir,
          job_log         => File::Spec->join($sample_obj->result_path, 'job-' . $step->name . '.yml'),
          pipeline        => $config->get('pipelines',                  $sample_obj->center),
          max_failed_runs => $config->get($project,                     'max_failed_runs'),
          out_dir         => $sample_obj->result_path,
          run_dir         => $run_dir,
          project_dir     => $project_dir,
          delay           => $delay,
          threads         => $procs,
          meta_id         => $job_meta->id,
          mapper_cmd      => File::Spec->join($project_dir, $PROGRAM_NAME),
          cluster         => $cluster,
          project         => $project,
        },
        gotcloud => {
          root    => $gotcloud_root,
          conf    => $gotcloud_conf,
          ref_dir => $gotcloud_ref,
          cmd     => File::Spec->join($gotcloud_root, 'gotcloud'),
        },
        sample => $sample_obj,
      },
      $job_file
      ) or die $tt->error();

    $logger->debug("wrote batch file to $job_file") if $debug;

    $jobs++;
    next if $self->app->global_options->{dry_run};

    my $job = CSG::Mapper::Job->new(cluster => $cluster);

    $logger->debug("submitting batch file $job_file") if $debug;

    try {
      $job->submit($job_file);
    }
    catch {
      if (not ref $_) {
        $logger->critical('Uncaught exception');
        $logger->debug($_) if $debug;

      } elsif ($_->isa('CSG::Mapper::Exceptions::Job::BatchFileNotFound')) {
        $logger->error($_->description);

      } elsif ($_->isa('CSG::Mapper::Exceptions::Job::BatchFileNotReadable')) {
        $logger->error($_->description);

      } elsif ($_->isa('CSG::Mapper::Exceptions::Job::SubmissionFailure')) {
        $logger->error($_->description);

      } elsif ($_->isa('CSG::Mapper::Exceptions::Job::ProcessOutput')) {
        $logger->error($_->description);
        $logger->debug($_->output) if $debug;

      } else {
        if ($_->isa('Exception::Class')) {
          chomp(my $error = $_->error);
          $logger->critical($error);
        } else {
          $logger->critical('something went sideways');
          print STDERR Dumper $_ if $debug;
        }
      }
    }
    finally {
      unless (@_) {
        $logger->info('submitted job (' . $job->job_id . ') for sample ' . $sample_obj->sample_id) if $verbose;

        $result->update({state_id => $schema->resultset('State')->find({name => 'submitted'})->id});
        $job_meta->update(
          {
            job_id       => $job->job_id(),
            submitted_at => $schema->now(),
          }
        );
      }
    };
  }
}

1;

__END__

=head1

CSG::Mapper::Command::launch - Launch remapping jobs
