package Test::CSG::Storage::SlotFS::Topmed;

use base qw(Test::CSG::Storage::Base);

use CSG::Base qw(file test formats);
use CSG::Storage::SlotFS::Topmed;

sub class {
  return 'CSG::Storage::SlotFS::Topmed';
}

sub startup : Test(startup => 2) {
  my ($self) = @_;

  my $fixtures = YAML::LoadFile(File::Spec->join($self->fixture_path, 'samples.yml'));

  for my $fixture (@{$fixtures}) {
    my $sample = $self->class->new(
      name     => $fixture->{sample_id},
      project  => $fixture->{project},
      prefix   => $fixture->{prefix},
    );

    isa_ok($sample, $self->class);
    push @{$self->{stash}->{samples}}, $sample;
  }
}

sub test_path : Test(2) {
  my ($self) = @_;
  for my $sample (@{$self->{stash}->{samples}}) {
    ok(-e $sample->path, 'slot path exists');
  }
}

sub test_initialize : Test(no_plan) {
  local $TODO = 'not implemented yet';
}

sub test_save_manifest : Test(8) {
  my ($self) = @_;

  for my $sample (@{$self->{stash}->{samples}}) {
    diag($sample->manifest);

    can_ok($sample, 'save_manifest');
    ok($sample->save_manifest, 'saved manifest');
    file_exists_ok($sample->manifest, 'manifest exists on disk');
    file_not_empty_ok($sample->manifest, 'manifest is not zero length');
  }
}

1;
