package CSG::Mapper::DB::Schema::ResultSet::Sample;

use base qw(DBIx::Class::ResultSet);
use CSG::Base;
use CSG::Mapper::Exceptions;

sub available_for {
  my ($self, $build, $step)  = @_;

  my $meth = qq{_available_for_$step};
  unless ($self->can($meth)) {
    CSG::Mapper::Exceptions::Sample::InvalidStep->throw('invalid availability');
  }

  return $self->$meth($build);
}

sub _available_for_bam2fastq {
  # TODO - needs to find samples where there is no result and where a result has
  #        is requested
  return shift->search(
    {
      -or => [
        {'results.id'  => undef},
        {'state.name' => 'requested'},
      ],
    }, {
      join     => {results => {results_states_steps => [qw(state step)]}},
      order_by => 'RAND()',
    }
  );
}

sub _available_for_align {
  my ($self, $build) = @_;

  my $schema     = $self->result_source->schema();
  my $candidates = $schema->resultset('ResultsStatesStep')->search(
    {
      'state.name' => 'completed',
      'step.name'  => 'bam2fastq',
    }, {
      join     => [qw(state step)],
      order_by => 'created_at desc',
    }
  );

  return map {$_->result->sample} $candidates->all();
}

1;
