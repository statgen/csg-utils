package Test::CSG::Storage::Base;

use base qw(Test::Class);

use CSG::Base qw(file test formats);
use CSG::Storage::Config;
use CSG::Storage::Slots::DB;

use Number::Bytes::Human qw(parse_bytes);

my $PREFIX = File::Temp::tempdir();

sub fixture_path {
  return qq{$FindBin::Bin/../t/fixtures};
}

sub prefix {
  return $PREFIX;
}

sub _startup : Test(startup => 1) {
  my ($self) = @_;
  diag('Using prefix ' . $self->prefix);
  ok(-e $self->prefix, 'prefix exists');
}

sub _shutdown : Test(shutdown => 1) {
  my ($self) = @_;
  diag('removing temporary directory ' . $self->prefix);
  remove_tree($self->prefix);
  ok(!-e $self->prefix, 'deleted prefix');
}

1;
