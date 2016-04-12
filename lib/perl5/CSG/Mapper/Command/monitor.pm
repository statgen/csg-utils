package CSG::Mapper::Command::monitor;

use CSG::Mapper -command;
use CSG::Base qw(file templates);
use CSG::Constants qw(:basic :mapping);
use CSG::Mapper::Config;
use CSG::Mapper::DB;

sub opt_spec {
  return (
    ['build=s', 'Reference build to run jobs against']
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

  unless ($opts->{build}) {
    $self->usage_error('Build is required');
  }

  unless ($opts->{build} =~ /37|38/) {
    $self->usage_error('Invalid build specified');
  }

  my $config = CSG::Mapper::Config->new(project => $self->app->global_options->{project});
  $self->{stash}->{config} = $config;
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $config = $self->{stash}->{config};
  my $logger = CSG::Mapper::Logger->new();

  my $cluster = $self->app->global_options->{cluster};
  my $project = $self->app->global_options->{project};
  my $debug   = $self->app->global_options->{debug};
  my $verbose = $self->app->global_options->{verbose};

  my $prefix      = $config->get($cluster, 'prefix');
  my $project_dir = qq{$FindBin::Bin/..};
  my $control_dir = File::Spec->join($project_dir, 'control');
  my $workdir     = $config->get($project, 'workdir');

  my $basedir = File::Spec->join($prefix, $workdir);
  $logger->debug("basedir: $basedir") if $debug;
  unless (-e $basedir) {
    make_path($basedir);
    $logger->debug('created basedir') if $debug;
  }

  my $log_dir = File::Spec->join($basedir, $config->get($project, 'log_dir'), 'monitor');
  $logger->debug("log_dir: $log_dir") if $debug;
  unless (-e $log_dir) {
    make_path($log_dir);
    $logger->debug('created log_dir') if $debug;
  }

  my $run_dir = File::Spec->join($basedir, $config->get($project, 'run_dir'));
  $logger->debug("run_dir: $run_dir") if $debug;
  unless (-e $run_dir) {
    make_path($run_dir);
    $logger->debug('created run_dir') if $debug;
  }

  my $job_file = File::Spec->join($run_dir, join($DASH, ('monitor', $cluster, $project, 'hg' . $opts->{build} . '.sh')));
  my $tt = Template->new(INCLUDE_PATH => qq($project_dir/templates/monitor));

  $tt->process(
    'all.sh.tt2', {
      settings => {
        cluster     => $cluster,
        project     => $project,
        project_dir => $project_dir,
        control_dir => $control_dir,
        mapper_cmd  => $PROGRAM_NAME,
        build       => $opts->{build},
      },
      job => {
        email   => $config->get($project, 'email'),
        account => $config->get($cluster, 'account'),
        workdir => $log_dir,
      },
    },
    $job_file
    )
    or croak $Template::ERROR;

  unless ($self->app->global_options->{dry_run}) {
    my $job = CSG::Mapper::Job->new(cluster => $cluster);

    try {
      $job->submit($job_file);
    }
    catch {
      if (not ref $_) {
        $logger->critical('Uncaught exception');
        $logger->debug($_) if $debug;

      } elsif ($_->isa('CSG::Mapper::Execption::Job::BatchFileNotFound')) {
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
        $logger->info('submitted job (' . $job->job_id . ") for monitoring project: $project on cluster: $cluster");
      }
    };
  }
}

1;

__END__

=head1

CSG::Mapper::Command::monitor - submit a monitor job to track samples and submit more jobs
