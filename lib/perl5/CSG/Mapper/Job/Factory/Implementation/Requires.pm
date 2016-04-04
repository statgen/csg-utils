package CSG::Mapper::Job::Factory::Implementation::Requires;

use Moose::Role;

requires(
  qw(
    job_id
    has_job_id
    elapsed
    state
    job_output_regexp
    job_submit_cmd
    job_kill_cmd
    )
);

1;
