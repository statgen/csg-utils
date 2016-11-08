package CSG::Mapper::Command::hold;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Mapper::DB;

my $schema = CSG::Mapper::DB->new();

sub opt_spec {
  return (
    ['sample=s', 'sample to work with', {required => 1}],
    ['step=s',   'what step is this sample currently processing', {required => 1}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $build   = $self->app->global_options->{build};
  my $cluster = $self->app->global_options->{cluster};

  my $step  = $schema->resultset('Step')->find({name => $opts->{step}});
  unless ($step) {
    $self->usage_error('Invalid step');
  }

  my $sample = $schema->resultset('Sample')->find({sample_id => $opts->{sample}});
  my $job    = $sample->jobs_for_build($build)->search({'step.name' => $step->name}, {join => 'step'})->first;
  if ($job->cluster ne $cluster) {
    say "The job that is processing sample $opts->{sample} is from a different cluster";
    exit 1;
  }

  $self->{stash}->{job} = $job;
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $meta    = $self->{stash}->{job};
  my $verbose = $self->app->global_options->{verbose};
  my $debug   = $self->app->global_options->{debug};

  my $logger  = CSG::Mapper::Logger->new(job_id => $meta->id);
  my $job     = CSG::Mapper::Job->new(
    cluster => $self->app->global_options->{cluster},
    job_id  => $opts->{job_id},
  );

  try {
    #$job->hold();
    $logger->info('holding job ' . $job->job_id);
  }
  catch {
    if (not ref $_) {
      $logger->critical('Uncaught exception');
      $logger->debug($_) if $debug;

    } elsif ($_->isa('CSG::Mapper::Exceptions::NoJobToCancel')) {
      $logger->error($_->description);

    } elsif ($_->isa('CSG::Mapper::Exceptions::Job::ProcessOutput')) {
      $logger->error($_->description);
      $logger->debug($_->output) if $debug;

    } elsif ($_->isa('CSG::Mapper::Exceptions::Job::CancellationFailure')) {
      $logger->error($_->description);

    }
  };
}

1;

__END__

=head1

CSG::Mapper::Command::hold - requeued and hold a running sample
