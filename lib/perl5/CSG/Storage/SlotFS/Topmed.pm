package CSG::Storage::SlotFS::Topmed;

use autodie qw(:all);
use Moose;

use Modern::Perl;
use Readonly;
use File::Path qw(make_path);
use IO::All;

our $VERSION = '0.1';

Readonly::Scalar my $PROJECT => 'topmed';
Readonly::Array my @PATHS    => (qw(incoming results logs run info));

has 'size' => (
  is      => 'rw',
  isa     => 'Int',
  default => sub {
    return 0;
  },
);

has 'project' => (
  is      => 'ro',
  isa     => 'Str',
  default => sub {
    return $PROJECT;
  },
);

for my $path (@PATHS) {
  has "${path}_path" => (
    is      => 'ro',
    isa     => 'Str',
    lazy    => 1,
    builder => "_build_${path}_path",
  );

  eval "sub _build_${path}_path { return File::Spec->join(shift->path, $path); }";
}

sub initialize {
  my ($self) = @_;

  my @skel_dirs = map {File::Spec->join($self->path, $_)} @PATHS;

  make_path(@skel_dirs, {error => \my $err});

  if (@{$err}) {
    my $errstr;
    for (@{$err}) {
      my ($key, $value) = %{$_};
      $errstr = $value if $key eq '';
    }

    CSG::Storage::Slots::Exceptions::Sample::FailedSkeletonDirectory->throw(error => $errstr);
  }

  return;
}

sub to_string {
  return shift->path;
}

sub manifest {
  File::Spec->join(shift->path, 'MANIFEST.yml');
}

sub save_manifest {
  my ($self) = @_;

  my $manifest = {slot => $self->slot_manifest};
  return YAML::Dump($manifest) > io($self->manifest);
}

with
  'CSG::Storage::SlotFS::Roles::Sample' => {-alias => {manifest => 'sample_manifest'}},
  'CSG::Storage::SlotFS::Roles::Slot'   => {-alias => {manifest => 'slot_manifest'}};

no Moose;
__PACKAGE__->meta->make_immutable;

1;
