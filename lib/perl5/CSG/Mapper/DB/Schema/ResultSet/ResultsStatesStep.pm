package CSG::Mapper::DB::Schema::ResultSet::ResultsStatesStep;

use base qw(DBIx::Class::ResultSet);
use CSG::Base;

sub current_results_by_step {
  my ($self, $build, $step) = @_;

  return $self->search(
    {
      'step.name'    => $step,
      'result.build' => $build,
    },
    {
      join     => [qw(result step)],
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search({}, {
      group_by => 'result_id'
  });
}

sub current_results_by_step_state {
  my ($self, $build, $step, $state) = @_;

  return $self->search(
    {},
    {
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search(
    {},
    {
      group_by => 'me.result_id',
    }
  )->as_subselect_rs->search(
    {
      'result.build' => $build,
      'state.name'   => $state,
      'step.name'    => $step,
    },
    {
      join => [qw(result step state)],
    }
  );
}

sub current_state_for_result {
  my ($self, $result_id) = @_;

  return $self->search(
    {
      'result.id' => $result_id,
    },
    {
      join     => 'result',
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search({}, {
      group_by => 'result_id'
  });
}

1;
