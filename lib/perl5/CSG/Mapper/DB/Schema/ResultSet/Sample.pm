package CSG::Mapper::DB::Schema::ResultSet::Sample;

use base qw(DBIx::Class::ResultSet);
use CSG::Base;

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

1;
