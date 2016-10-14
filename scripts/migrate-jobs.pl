#!/usr/bin/env perl

use CSG::Base;
use CSG::Mapper::DB;

my $schema = CSG::Mapper::DB->new();

for my $job ($schema->resultset('Job')->all()) {
  $schema->resultset('ResultsStatesStep')->search({result_id => $job->result_id})->update({job_id => $job->id});
}
