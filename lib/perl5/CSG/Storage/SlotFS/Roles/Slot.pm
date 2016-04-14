package CSG::Storage::SlotFS::Roles::Slot;

use autodie qw(:all);
use Moose::Role;

use Modern::Perl;
use File::Spec;

use CSG::Storage::Slots;
use CSG::Storage::Types;

requires qw(size to_string project);

has 'name' => (
  is       => 'ro',
  isa      => 'Str',
  required => 1,
);

has 'prefix' => (
  is      => 'ro',
  isa     => 'ValidPrefixPath',
  default => sub {
    return '/net'
  },
);

has 'path' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_path',
);

has 'slot' => (
  is      => 'ro',
  isa     => 'CSG::Storage::Slots',
  lazy    => 1,
  builder => '_build_slot',
);

sub _build_slot {
  my ($self) = @_;

  return CSG::Storage::Slots->find_or_create(
    name    => $self->name,
    size    => $self->size,
    project => $self->project,
  );
}

sub _build_path {
  my ($self) = @_;
  return File::Spec->join($self->prefix, $self->slot->path);
}

1;
