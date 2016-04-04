package CSG::Mapper::Job::Factory::Implementation::csg;

use CSG::Base qw(cmd file);
use CSG::Constants;
use CSG::Mapper::Util qw(parse_time);

use Moose;

Readonly::Scalar my $PREFIX     => q{/usr/cluster/bin};
Readonly::Scalar my $SACCT_CMD  => File::Spec->join($PREFIX, 'sacct');
Readonly::Scalar my $SQUEUE_CMD => File::Spec->join($PREFIX, 'squeue');
Readonly::Scalar my $SBATCH_CMD => File::Spec->join($PREFIX, 'sbatch');
Readonly::Scalar my $SCANEL_CMD => File::Spec->join($PREFIX, 'scancel');

Readonly::Scalar my $JOB_ELAPSED_TIME_FORMAT   => $SACCT_CMD . q{ -j %d -X -n -o elapsed};
Readonly::Scalar my $JOB_STATE_CMD_FORMAT      => $SACCT_CMD . q{ -j %d -X -n -o state%%20};
Readonly::Scalar my $JOB_TIME_REMAINING_FORMAT => $SQUEUE_CMD . q{ -h -o %%L -j %d};
Readonly::Scalar my $JOB_OUTPUT_REGEXP         => qr/^Submitted batch job (?<jobid>\d+)$/i;

Readonly::Hash my %JOB_STATES => (
  RUNNING   => 'running',
  COMPLETED => 'completed',
  FAILED    => 'failed',
  REQUEUED  => 'requeued',
  CANCELLED => 'cancelled',
);

has 'job_id'            => (is => 'rw', isa => 'Int',       predicate => 'has_job_id');
has 'job_output_regexp' => (is => 'ro', isa => 'RegexpRef', default   => sub {return $JOB_OUTPUT_REGEXP});
has 'job_submit_cmd'    => (is => 'ro', isa => 'Str',       default   => sub {return $SBATCH_CMD});
has 'job_submit_kill'   => (is => 'ro', isa => 'Str',       default   => sub {return $SCANCEL_CMD});

sub elapsed {
  my ($self) = @_;
  my $cmd = sprintf $JOB_ELAPSED_TIME_FORMAT, $self->job_id;
  chomp(my $time = capture(EXIT_ANY, $cmd));
  return parse_time($time);
}

sub state {
  my ($self) = @_;
  my $cmd = sprintf $JOB_STATE_CMD_FORMAT, $self->job_id;
  chomp(my $state = capture(EXIT_ANY, $cmd));
  $state =~ s/^\s+|\s+$//g;
  return $JOB_STATES{$state};
}

sub _time_remaining {
  return capture(sprintf($JOB_TIME_REMAINING_FORMAT, shift->job_id));
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

