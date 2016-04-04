package CSG::Mapper::Command::resubmit;

use CSG::Mapper -command;
use CSG::Base qw(file);
use CSG::Constants qw(:mapping);
use CSG::Mapper::DB;
use CSG::Mapper::Logger;
use CSG::Mapper::Sample;

sub opt_spec {
  return (
    ['build|b=i',     'Reference build to use (valid values: 37 or 38)'],
    ['sample-id|s=s', 'Sample ID to resubmit (e.g. NWD12345)'],
    ['failed|f',      'Resubmit all failed jobs'],
    ['do-not-purge',  'Prevent the deletion of OUT_DIR'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($self->app->global_options->{cluster}) {
    $self->usage_error('cluster is required');
  }

  unless ($self->app->global_options->{cluster} =~ /$VALID_CLUSTER_REGEXPS/) {
    $self->usage_error('Invalid cluster environment');
  }

  unless ($opts->{build}) {
    $self->usage_error('reference build is required');
  }

  unless ($opts->{build} =~ /37|38/) {
    $self->usage_error('invalid reference build');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my @samples = ();
  my $logger  = CSG::Mapper::Logger->new();
  my $schema  = CSG::Mapper::DB->new();
  my $state   = $schema->resultset('State')->find({name => 'requested'});
  my $debug   = $self->app->global_options->{debug};

  if ($opts->{sample_id}) {
    push @samples, $opts->{sample_id};
  } elsif ($opts->{failed}) {
    my $results = $schema->resultset('Result')->search(
      {
        'me.build' => $opts->{build},
        'state_id' => $schema->resultset('State')->find({name => 'failed'})->id,
      }
    );

    @samples = map {$_->sample->sample_id} $results->all();
  }

  for my $sample (@samples) {
    my $result = $schema->resultset('Result')->search(
      {
        'me.build'         => $opts->{build},
        'sample.sample_id' => $sample,
      }, {
        join => 'sample',
      }
    );

    if ($result->count > 1) {
      croak 'Found multiple results for this sample/build';
    }

    $result->first->update(
      {
        state_id => $state->id,
      }
    );

    $logger->info("marked sample[$sample] for remapping");

    my $sample_obj = CSG::Mapper::Sample->new(
      cluster => $self->app->global_options->{cluster},
      record  => $result->first->sample,
      build   => $opts->{build},
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
