#!/usr/bin/env perl

use CSG::Base qw(file);
use CSG::Mapper::DB;

my @tmp_dirs = (
  qw(
    /net/topmed6/working/mapping/tmp/topmed/hg38
    /net/topmed7/working/mapping/tmp/topmed/hg38
    /net/topmed8/working/mapping/tmp/topmed/hg38
    )
);

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
  my $fastq = $job->results_states_steps->first->result->sample->fastqs->first;
  my $path  = undef;

  if ($fastq) {
    $path = dirname($fastq->path);
  } else {
    my $sample_id = $job->results_states_steps->first->result->sample->sample_id;
    for (@tmp_dirs) {
      my $tmp = File::Spec->join($_, $sample_id);
      if (-e $tmp) {
        $path = $tmp;
        last;
      }
    }
  }

  say $job->id . ' => ' . $job->tmp_dir . ' => ' . $path;
  $job->update({tmp_dir => $path});
}
