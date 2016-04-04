package CSG::Mapper::Job;

use Moose;

use CSG::Base qw(cmd);
use CSG::Constants qw(:mapping);
use CSG::Types;
use CSG::Mapper::Exceptions;
use CSG::Mapper::Job::Factory;
use CSG::Mapper::Logger;

has 'cluster' => (is => 'ro', isa => 'ValidCluster', required  => 1);
has 'job_id'  => (is => 'rw', isa => 'Int',          predicate => 'has_job_id');
has 'factory' => (
  is      => 'ro',
  isa     => 'ValidJobFactory',
  lazy    => 1,
  builder => '_build_factory',
  handles => [
    qw(
      elapsed
      state
      job_output_regexp
      job_submit_cmd
      job_kill_cmd
      _time_remaining
      )
  ],
);

sub _build_factory {
  my ($self) = @_;
  my $class  = __PACKAGE__ . q{::Factory};
  my $opts   = {};

  if ($self->has_job_id) {
    $opts->{job_id} = $self->job_id;
  }

  return $class->create($self->cluster, $opts);
}

sub elapsed_seconds {
  my $e = shift->elapsed;
  return ($e->days * 24 * 3600) + ($e->hours * 3600) + ($e->minutes * 60) + $e->seconds;
}

sub submit {
  my ($self, $file) = @_;

  CSG::Mapper::Exceptions::Job::BatchFileNotFound->throw() unless -e $file;
  CSG::Mapper::Exceptions::Job::BatchFileNotReadable->throw() unless -r $file;

  try {
    chomp(my $output = capture($self->job_submit_cmd, $file));
    my $regexp = $self->job_output_regexp;
    if ($output =~ /$regexp/) {
      $self->job_id($+{jobid});
    } else {
      CSG::Mapper::Exceptions::Job::ProcessOutput->throw(output => $output);
    }
  } catch {
    die $_ unless blessed $_ and $_->can('rethrow');

    $_->rethrow if $_->can('rethrow');

    CSG::Mapper::Exceptions::Job::SubmissionFailure->throw(error => $_);
  };

  return;
}

sub cancel {
  my ($self) = @_;

  CSG::Mapper::Exceptions::Job::NoJobToCancel->throw() unless $self->has_job_id;

  try {
    chomp(my $output = capture($self->job_kill_cmd, $self->job_id));

    if ($output) {
      CSG::Mapper::Exceptions::Job::ProcessOutput->throw(output => $output);
    }
  } catch {
    die $_ unless blessed $_ and $_->can('rethrow');

    $_->rethrow if $_->can('rethrow');

    CSG::Mapper::Exceptions::Job::CancellationFailure->throw(error => $_);
  };

  return;
}

sub time_remaining {
  my ($self) = @_;

  for my $regexp (@TIME_FORMAT_REGEXPS) {
    if ($self->_time_remaining() =~ $regexp) {
      return (($+{days} * 24) + $+{hours}) if $+{days} and $+{hours};
      return $+{hours} if $+{hours};
      return int($+{seconds} / 60 / 60) if $+{seconds};
    }
  }

  return 0;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
