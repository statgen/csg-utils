#!/usr/bin/env perl

use FindBin qw($Bin);
use lib qq($Bin/../lib/perl5);
use DBIx::Class::Schema::Loader qw(make_schema_at);

use CSG::Storage::Config;

my $conf = CSG::Storage::Config->new();

make_schema_at(
  'CSG::Storage::Slots::DB::Schema', {
    debug          => 1,
    dump_directory => qq($Bin/../lib/perl5),
    components     => [qw(InflateColumn::DateTime CSG::CreatedAt)],
  },
  [$conf->dsn, $conf->db_user, $conf->db_pass]
);
