package CSG::Mapper::DB::Schema::ResultSet::Result;

use base qw(DBIx::Class::ResultSet);
use CSG::Base;

sub completed {return shift->_state('completed', @_);}
sub failed    {return shift->_state('failed',    @_);}
sub cancelled {return shift->_state('cancelled', @_);}
sub requested {return shift->_state('requested', @_);}
sub started   {return shift->_state('started',   @_);}
sub submitted {return shift->_state('submitted', @_);}

sub _state {
  my ($self, $state, $step, $build) = @_;

  return $self->search(
    {
      'me.build'   => $build,
      'state.name' => $state,
      'step.name'  => $step,
    }, {
      join     => {results_states_steps => [qw(state step)]},
      group_by => 'me.sample_id',
    }
  );
}

1;
