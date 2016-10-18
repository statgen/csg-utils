#!/usr/bin/env perl

use CSG::Base qw(file);
use CSG::Mapper::DB;

my $schema = CSG::Mapper::DB->new();
my $jobs   = $schema->resultset('Job')->search(
  {
    tmp_dir => {like => '/tmp/%'},
    job_id  => {'!=' => 0},
  }, {
    join     => 'results_states_steps',
    group_by => 'job_id',
  }
);

for my $job ($jobs->all) {
  my $path = dirname($job->results_states_steps->first->result->sample->fastqs->first->path);
  say $job->id . ' => ' . $job->tmp_dir. ' => ' . $path;
  $job->update({tmp_dir => $path});
}
