#!/usr/bin/env perl

use CSG::Base qw(file);
use CSG::Mapper::DB;

my $schema = CSG::Mapper::DB->new();

for my $job ($schema->resultset('Job')->all) {
  my $result = $job->result;
  next unless $result;

  my $sample = $job->result->sample;
  next unless $sample->fastqs->count;
  next unless $job->results_states_steps->count;

  my $tmp_dir = dirname($sample->fastqs->first->path);

  $job->update({tmp_dir => $tmp_dir});
  say 'set job ' . $job->id . " tmp_dir to $tmp_dir";
}
