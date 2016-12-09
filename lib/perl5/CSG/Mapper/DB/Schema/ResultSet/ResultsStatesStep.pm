package CSG::Mapper::DB::Schema::ResultSet::ResultsStatesStep;

use base qw(DBIx::Class::ResultSet);
use CSG::Base;

sub current_results_by_step {
  my ($self, $build, $step) = @_;

  return $self->search(
    {
      'step.name'    => $step,
      'result.build' => $build,
      'project.name' => 'topmed',
    },
    {
      join     => ['step', {result => {sample => 'project'}}],
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search({}, {
      group_by => 'result_id'
  });
}

sub current_results_by_step_state {
  my ($self, $build, $step, $state) = @_;

  return $self->search(
    {
      'result.build' => $build,
      'project.name' => 'topmed',
      'step.name'    => $step,
    },
    {
      join     => ['step', {result => {sample => 'project'}}],
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search(
    {},
    {
      group_by => 'me.result_id',
    }
  )->as_subselect_rs->search(
    {
      'state.name'   => $state,
    },
    {
      join => [qw(state)],
    }
  );
}

sub current_state_for_result {
  my ($self, $result_id) = @_;

  return $self->search(
    {
      'result.id'    => $result_id,
      'project.name' => 'topmed',
    },
    {
      join     => {result => {sample => 'project'}},
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search({}, {
      group_by => 'result_id'
  });
}

sub running_by_step {
  my ($self, $build, $step) = @_;

  return $self->search(
    {
      'step.name'    => $step,
      'result.build' => $build,
      'state.name'   => {'=' => [qw(submitted started)]},
      'project.name' => 'topmed',
    },
    {
      join     => [qw(step state), {result => {sample => 'project'}}],
      order_by => 'created_at desc',
    }
  )->as_subselect_rs->search({}, {
      group_by => 'result_id'
  });
}

1;
