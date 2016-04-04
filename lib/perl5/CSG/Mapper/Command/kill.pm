package CSG::Mapper::Command::kill;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Mapper::DB;

sub opt_spec {
  return (['job-id=s', 'job to kill based on the clusters assigned id'],);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($self->app->global_options->{cluster}) {
    $self->usage_error('Cluster is required');
  }

  unless ($opts->{job_id}) {
    $self->usage_error('Job ID is required');
  }

  my $schema = CSG::Mapper::DB->new();
  my $job = $schema->resultset('Job')->find({job_id => $opts->{job_id}});
  unless ($job) {
    $self->usage_error("Job ID, $opts->{job_id}, does not exist");
  }

  if ($job->cluster ne $self->app->global_options->{cluster}) {
    $self->usage_error("Job ID, $opts->{job_id}, does not belong to specified cluster");
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
    $job->cancel();
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
  }
  finally {
    unless (@_) {
      $logger->info('cancelled job ' . $job->job_id) if $verbose;
      $meta->cancel();
    }
  };
}

1;

__END__

=head1

CSG::Mapper::Command::kill - cancel a remapping job
