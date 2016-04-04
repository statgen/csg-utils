#!/usr/bin/env perl

use FindBin qw($Bin);
use local::lib qq($Bin/../local);
use lib qq($Bin/../lib/perl5);

use CSG::Base qw(parsers file);
use CSG::Constants qw(:basic :mapping);
use CSG::Mapper::Config;
use CSG::Mapper::DB;

my $schema = CSG::Mapper::DB->new();
my $csv    = Class::CSV->parse(
  filehandle => io->stdin->tie,
  fields     => [qw(center run_dir filename study pi sample_id state_b37 state_b38 fullpath)],
);

my @lines = @{$csv->lines()};
shift @lines;

for my $line (@lines) {
  next if $line->state_b37 != 20;

  my $state = $schema->resultset('State')->find({name => 'completed'});
  my $sample = $schema->resultset('Sample')->search({sample_id => $line->sample_id})->first;

  next unless $sample;

  my $result = $schema->resultset('Result')->find_or_create(
    {
      sample_id => $sample->id,
      state_id  => $state->id,
      build     => '37',
    }
  );
}
