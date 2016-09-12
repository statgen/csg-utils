## no critic (NamingConventions::Capitalization, Subroutines::RequireFinalReturn)
package CSG::Mapper::Command::update;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Constants qw(:basic :mapping);
use CSG::Mapper::DB;
use CSG::Mapper::Exceptions;
use CSG::Mapper::Logger;

my $schema = CSG::Mapper::DB->new();

sub opt_spec {
  return (
    ['meta-id=i',    'Job meta data db record id', {required => 1}],
    ['start',        'Mark a sample started'],
    ['job-id=i',     'Add the clusters job id for a given sample'],
    ['node=s',       'Update what node(s) a sample is running on in the cluster'],
    ['state=s',      'Update the jobs state (valid states: failed|submitted|completed|cancelled|requested)'],
    ['exit-code=i',  'Update the exit code from a given sample'],
    ['step=s',       'Job step [bam2fastq|align|cloud-align|all]'],
    ['fastq-list=s', 'List of fastq files generated by bamUtil'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $meta = $schema->resultset('Job')->find($opts->{meta_id});
  unless ($meta) {
    $self->usage_error('unable to locate the job meta data record');
  }

  if ($opts->{state}) {
    my $state = $schema->resultset('State')->find({name => $opts->{state}});
    unless ($state) {
      $self->usage_error('invalid job state');
    }

    $self->{stash}->{state} = $state;
  }

  if ($opts->{start} and $meta->started_at) {
    $self->usage_error('job has already started');
  }

  if (defined $opts->{exit_code} and $meta->ended_at) {
    $self->usage_error('job has already ended');
  }

  if ($opts->{step}) {
    my $step = $schema->resultset('Step')->find({name => $opts->{step}});
    unless ($step) {
      $self->usage_error('invalid job step');
    }

    $self->{stash}->{step} = $step;
  }

  if ($opts->{fastq_list}) {
    unless (-e $opts->{fastq_list}) {
      $self->usage_error('fastq list does not exist');
    }
  }

  $self->{stash}->{meta} = $meta;
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $meta   = $self->{stash}->{meta};
  my $state  = $self->{stash}->{state};
  my $step   = $self->{stash}->{step};
  my $logger = CSG::Mapper::Logger->new(job_id => $meta->id);
  my $params = {};

  if ($opts->{start}) {
    $params->{started_at} = $schema->now();
  }

  if ($opts->{job_id}) {
    $params->{job_id} = $opts->{job_id};
  }

  if ($opts->{node}) {
    $params->{node} = $opts->{node};
  }

  if ($opts->{step}) {
    $params->{step_id} = $step->id;
  }

  if (defined $opts->{exit_code}) {
    $params->{exit_code} = $opts->{exit_code};
    $params->{ended_at}  = $schema->now();
  }

  if (keys %{$params}) {
    for my $key (keys %{$params}) {
      $logger->info("updating $key to $params->{$key}");
    }

    $meta->update($params);
  }

  if ($state) {
    $logger->info('changing result state from ' . $meta->result->state->name . ' to ' . $state->name);
    $meta->result->update({state_id => $state->id});
  }

  if ($opts->{fastq_list}) {
    for my $line (io($opts->{fastq_list})->chomp->getlines) {
      next if $line =~ /^MERGE/;
      my ($sample_id, $fastq1, $fastq2, @rg) = split(/$TAB/, $line);

      next if $meta->result->sample->has_fastq($fastq1);

      unless ($sample_id eq $meta->result->sample->sample_id) {
        CSG::Mapper::Exceptions::Sample::FastqMismatch->throw();
      }

      unless (-e $fastq1) {
        CSG::Mapper::Exceptions::Sample::FastqNotFound->throw();
      }

      unless ($rg[0] eq '@RG') {
        CSG::Mapper::Exceptions::Sample::Fastq::MissingRG->throw();
      }

      # XXX - It appears read group headers are not the same across
      #       samples making this test blow up when it shouldn't. might
      #       bring it back if I can pin down these headers more.
      #
      # my $headers = {map {(split(/$COLON/, $_))[0] => 1} @rg[1 .. $#rg]};
      # for (@READ_GROUP_FIELDS) {
      #   unless (exists $headers->{$_}) {
      #     CSG::Mapper::Exceptions::Sample::Fastq::MissingHeader->throw(header => $_);
      #   }
      # }

      $meta->result->sample->add_to_fastqs({
        build      => $meta->result->build,
        path       => $fastq1,
        read_group => join('\t', @rg),
      });

      $logger->info('added fastq ' . $fastq1 . ' to sample ' . $meta->result->sample->sample_id);
    }
  }
}

1;

__END__

=head1

CSG::Mapper::Command::update - update remapping jobs
