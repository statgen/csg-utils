package CSG::Mapper::DB::Schema::ResultSet::ResultsStatesStep;

use base qw(DBIx::Class::ResultSet);
use CSG::Base;

sub current_results_by_step {
  my ($self, $step) = @_;

  return $self->search(
    {
      'step.name' => $step,
    },
    {
      join     => 'step',
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search({}, {
      group_by => 'result_id'
  });
}

sub current_results_by_step_state {
  my ($self, $step, $state) = @_;

  return $self->search(
    {
      'step.name'  => $step,
      'state.name' => $state,
    },
    {
      join     => [qw(step state)],
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search({}, {
      group_by => 'result_id'
  });
}

1;
