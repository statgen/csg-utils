#!/usr/bin/env perl

use FindBin qw($Bin);
use DBIx::Class::Schema::Loader qw(make_schema_at);
use CSG::Base;
use CSG::Mapper::Config;

my $config = CSG::Mapper::Config->new(project => 'topmed');

make_schema_at(
  'CSG::Mapper::DB::Schema', {
    debug          => 1,
    dump_directory => qq($Bin/../lib/perl5),
    components     => [qw(InflateColumn::DateTime CSG::CreatedAt)],
  },
  [$config->dsn, $config->get('db','user'), $config->get('db','pass')]
);
