package CSG::Mapper::Command::resubmit;


use CSG::Mapper -command;
use CSG::Base qw(file);
use CSG::Constants qw(:mapping);
use CSG::Mapper::DB;
use CSG::Mapper::Logger;
use CSG::Mapper::Sample;

sub opt_spec {
  return (
    ['sample-id|s=s', 'Sample ID to resubmit (e.g. NWD12345)'],
    ['failed|f',      'Resubmit all failed jobs'],
    ['do-not-purge',  'Prevent the deletion of OUT_DIR'],
    ['step=s',        'Job step to resubmit (e.g. bam2fastq, cloud-align)', {required => 1}],
  );
}

sub execute {
  my ($self, $opts, $args) = @_;

  my @samples = ();
  my $logger  = CSG::Mapper::Logger->new();
  my $schema  = CSG::Mapper::DB->new();
  my $debug   = $self->app->global_options->{debug};
  my $build   = $self->app->global_options->{build};
  my $state   = $schema->resultset('State')->find({name => 'requested'});
  my $step    = $schema->resultset('Step')->find({name => $opts->{step}});

  unless ($step) {
    $logger->critical("Step, $opts->{step}, is not valid");
    exit 1;
  }

  if ($opts->{sample_id}) {
    push @samples, $schema->resultset('Sample')->find({sample_id => $opts->{sample_id}})->id;
  } elsif ($opts->{failed}) {
    my $results = $schema->resultset('ResultsStatesStep')->current_results_by_step_state($build, $opts->{step}, 'failed');
    @samples = map {$_->result->sample->id} $results->all();
  }

  for my $sample (@samples) {
    my $result = $schema->resultset('Result')->find({build => $build, sample_id => $sample});

    unless ($result) {
      $logger->critical("No result found for sample[$sample] in build[$build]");
      exit 1;
    }

    $result->add_to_results_states_steps({
        state_id => $state->id,
        step_id  => $step->id,
    });

    $logger->info("marked sample[$sample] for remapping");

    my $sample_obj = CSG::Mapper::Sample->new(
      cluster => $self->app->global_options->{cluster},
      record  => $result->sample,
      build   => $build,
    );

    $logger->debug("SAMPLE[$sample] OUT_DIR: " . $sample_obj->result_path) if $debug;

    unless ($opts->{do_not_purge}) {
      if (-e $sample_obj->result_path) {
        $logger->info("removing existing OUT_DIR for sample[$sample]");
        remove_tree($sample_obj->result_path);
      }
    }
  }
}

1;

__END__

=head1

CSG::Mapper::Command::resubmit - resubmit remapping jobs
