package CSG::Storage::SlotFS::Roles::Slot;

use autodie qw(:all);
use Moose::Role;

use Modern::Perl;
use File::Spec;

use CSG::Storage::Slots;
use CSG::Types;
use CSG::Base qw(file);

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
    return '/net';
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

has 'exclude' => (
  is      => 'ro',
  isa     => 'Maybe[Int]',
  default => sub {
    return undef;
  },
);

around 'path' => sub {
  my ($orig, $self) = @_;

  unless (-e $self->$orig) {
    make_path($self->$orig);
  }

  return $self->$orig;
};

sub _build_slot {
  my ($self) = @_;

  return CSG::Storage::Slots->find_or_create(
    name    => $self->name,
    size    => $self->size,
    project => $self->project,
    exclude => $self->exclude,
  );
}

sub _build_path {
  my ($self) = @_;
  return File::Spec->join($self->prefix, $self->slot->path);
}

sub find {
  my ($self, %params) = @_;

  return unless exists $params{name} and exists $params{project};

  return CSG::Storage::Slots->find(
    name    => $params{name},
    project => $params{project},
  );
}

sub manifest {
  my ($self) = @_;

  my $pool    = $self->slot->_record->pool;
  my $type    = $pool->type;
  my $project = $pool->project;

  return {
    name    => $self->name,
    sha1    => $self->slot->sha1,
    size    => $self->size,
    pool    => {
      id       => $pool->id,
      name     => $pool->name,
      hostname => $pool->hostname,
      project  => {
        id   => $project->id,
        name => $project->name,
      },
      type => {
        id   => $type->id,
        name => $type->name,
      }
    },
  };
}

1;
