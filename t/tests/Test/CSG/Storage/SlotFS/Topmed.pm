package Test::CSG::Storage::SlotFS::Topmed;

use base qw(Test::CSG::Storage::Base);
use Test::More;

use Modern::Perl;
use File::Stat;
use YAML qw(LoadFile);

use CSG::Storage::SlotFS::Topmed;

sub class {
  return 'CSG::Storage::SlotFS::Topmed';
}

sub startup : Test(startup => 2) {
  my ($self) = @_;

  my $fixtures = LoadFile(File::Spec->join($self->fixture_path, 'samples.yml'));

  for my $fixture (@{$fixtures}) {
    my $sample = $self->class->new(
      name     => $fixture->{sample_id},
      project  => $fixture->{project},
      prefix   => $fixture->{prefix},
    );

    isa_ok($sample, $self->class);
    diag $sample->path;

    push @{$self->{stash}->{samples}}, $sample;
  }
}

sub test_path : Test(no_plan) {
  local $TODO = 'not implemented yet';
}

sub test_initialize : Test(no_plan) {
  local $TODO = 'not implemented yet';
}

1;
