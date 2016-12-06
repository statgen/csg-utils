package CSG::Mapper::Command::align;

use CSG::Mapper -command;
use CSG::Base qw(cmd);
use CSG::Mapper::Config;
use CSG::Mapper::DB;
use CSG::Mapper::Logger;

use Parallel::ForkManager;

my $schema = CSG::Mapper::DB->new();
my $logger = CSG::Mapper::Logger->new();

sub opt_spec {
  return (
    ['limit|l=i',    'total number of samples to process',    {required => 1}],
    ['concurrent=i', 'number of concurrent samples to align', {required => 1}],
    ['requested',    'rerun samples in requested state'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($self->app->global_options->{cluster} ne 'dummy') {
    $self->usage_error('only cluster type "dummy" is accepted');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $jobs       = 0;
  my $parent_pid = $PROCESS_ID;

  my $cluster    = $self->app->global_options->{cluster};
  my $project    = $self->app->global_options->{project};
  my $build      = $self->app->global_options->{build};
  my $debug      = $self->app->global_options->{debug};
  my $verbose    = $self->app->global_options->{verbose};

  my $config     = CSG::Mapper::Config->new(project => $project);
  my $pm         = Parallel::ForkManager->new($opts->{concurrent});
  my $step       = $schema->resultset('Step')->find({name => 'cloud-align'});

  # XXX - LOCK TABLES 'samples';
  my @samples = ();
  my $search = {};
  my $attrs  = {
    order_by => 'RAND()',
    for      => 'update',
  };

  try {
    sleep(10);
    $schema->txn_do(
      sub {
        my $sample_rs = $schema->resultset('Sample')->search($search, $attrs);

        for my $sample ($sample_rs->all) {
          last if $opts->{limit} and $jobs >= $opts->{limit};

          unless ($sample->is_available($step->name, $build)) {
            $logger->info('sample ' . $sample->sample_id . ' is not available for processing') if $debug;
            next;
          }

          if ($opts->{requested}) {
            unless ($sample->is_requested($step->name, $build)) {
              $logger->info('sample ' . $sample->sample_id . ' is not in requested state') if $debug;
              next;
            }
          } else {
            if ($sample->is_requested($step->name, $build)) {
              $logger->info('sample ' . $sample->sample_id . ' is in requested state, skipping') if $debug;
              next;
            }
          }

          unless ($sample->fastqs->count) {
            $logger->debug('no fastq files recorded for sample ' . $sample->sample_id) if $debug;
            next;
          }

          if ($sample->has_fastqs_with_unpaired_reads) {
            $logger->debug('found unpaired reads in the fastqs for sample ' . $sample->sample_id) if $debug;
            next;
          }

          my $result = $sample->result_for_build($build);

          $schema->txn_do(
            sub {

              my $sample_ts = $schema->resultset('Sample')->search(
                {
                  id => $sample->id,
                },
                {
                  for => 'update',
                }
              );

              die if $sample_ts->count > 1;

              my $result_ts = $sample_ts->first->result_for_build($build);
              $result_ts->add_to_results_states_steps(
                {
                  state_id => $schema->resultset('State')->find({name => 'requested'})->id,
                  step_id  => $step->id,
                }
              );
            }
          );

          push @samples, $sample;

          $jobs++;
        }
      }
    );
  } catch {
    print Dumper $_;
    $logger->critical('transaction failed');
  };
  # XXX _ UNLOCK TABLES 'samples';

=cut
  exit;
  for my $sample (@samples) {
    my $sample_obj = CSG::Mapper::Sample->new(
      cluster => $cluster,
      record  => $sample,
      build   => $build,
    );

    my $child_pid = $pm->start and next;

    my $run_dir   = $sample_obj->state_dir;
    my $out       = sprintf 'cloud-align-%d_%d.out', $parent_pid, $child_pid;
    my $err       = sprintf 'cloud-align-%d_%d.err', $parent_pid, $child_pid;

    my $out_fh = IO::File->new(File::Spec->join($run_dir, $out), 'w+');
    my $err_fh = IO::File->new(File::Spec->join($run_dir, $err), 'w+');

    capture {
      my $result = $sample->result_for_build;
      my $status = $result->current_status;
      my $job    = $status->job; # FIXME - no job record yet

      $job->update(
        {
          job_id => $child_pid,
        }
      );

      next unless $status->state->name 'started';

      $result->add_to_results_states_steps(
        {
          state_id => $started_id,
          step_id  => $step_id,
          job_id   => $job->id, # FIXME - not sure what job record to use here
        }
      );

      # pid = job_id-node-ctime
      # my $script = `mapper launch --step cloud-align --sample $sample->sample_id --cluster dummy --requested`;
      #
      # XXX
      #   - job_id
      #   - sample_id
      # my $cmd = "sh $script";
      $ENV{MAPPER_JOB_ID} = $child_pid;
      run(EXIT_ANY, $cmd);
      # TODO - run gce-align.sh script or maybe a generated
      #        or maybe a generated cloud-align.sh batch script.

    } stdout => $out_fh, stderr => $err_fh;

    $pm->finish;
  }

  $pm->wait_all_children;
=cut
}

1;

__END__

=head1

CSG::Mapper::Command::align - run alignment jobs
