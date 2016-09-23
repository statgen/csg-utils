package CSG::Mapper::DB::Schema::ResultSet::Sample;

use base qw(DBIx::Class::ResultSet);
use CSG::Base;
use CSG::Mapper::Exceptions;

sub ready_for_alignment {
  my ($self, $build) = @_;

  return $self->search(
    {
      'step.name'           => 'bam2fastq',
      'state.name'          => 'submitted',
      'results.build'       => $build,
      'results.exported_at' => undef,
      'jobs.submitted_at'   => {'!=' => undef},
      'jobs.started_at'     => {'!=' => undef},
      'jobs.ended_at'       => {'!=' => undef},
      'jobs.exit_code'      => 0,
    }, {
      join => {'results' => ['state', {'jobs' => 'step'}]},
      '+select' => ['results.id', 'jobs.id'],
      '+as'     => [qw(result_id job_id)],
    }
  );
}

sub available_for {
  my ($self, $build, $step)  = @_;

  # TODO - this still isn't quite right in some cases
  #
  my $meth = qq{_available_for_$step};
  unless ($self->can($meth)) {
    CSG::Mapper::Exceptions::Sample::InvalidStep->throw('invalid availability');
  }

  return $self->$meth($build);
}

sub _available_for_bam2fastq {
  return shift->search(
    {
      'results.id' => undef,
    }, {
      join     => 'results',
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
