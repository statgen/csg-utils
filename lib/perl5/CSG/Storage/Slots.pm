package CSG::Storage::Slots;

use CSG::Base;
use CSG::Constants;
use CSG::Storage::Slots::DB;
use CSG::Storage::Slots::Exceptions;
use CSG::Types;

use Moose;
use Digest::SHA qw(sha1_hex);
use overload '""' => sub {shift->to_string};
use Path::Tiny ();

has '_record' => (
  is        => 'rw',
  isa       => __PACKAGE__ . '::DB::Schema::Result::Slot',
  predicate => 'has_record',
  handles   => [
    qw(
      sha1
      pool_id
      )
  ],
);

has 'name' => (
  is       => 'ro',
  isa      => 'Str',
  required => 1,
  trigger  => \&_set_name
);

has 'project' => (
  is       => 'ro',
  isa      => 'Str',
  required => 1,
  trigger  => \&_set_project
);

has 'size' => (
  is       => 'rw',
  isa      => 'ValidSlotSize',
  required => 1,
  trigger  => \&_set_size
);

has 'path' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_path'
);

has 'prefix' => (
  is      => 'ro',
  isa     => 'ValidPrefixPath',
  default => sub {
    return '/net';
  },
);

has 'parent' => (
  is        => 'rw',
  isa       => 'Maybe[Str]',
  predicate => 'has_parent',
  default   => sub {
    return undef;
  }
);

around [qw(name project size path sha1)] => sub {
  my $orig = shift;
  my $self = shift;

  unless ($self->has_record) {
    my $schema = CSG::Storage::Slots::DB->new();

    my $project = $schema->resultset('Project')->find({name => $self->{project}});
    unless ($project) {
      CSG::Storage::Slots::Exceptions::Project::DoesNotExist->throw();
    }

    if (__PACKAGE__->exists(name => $self->{name}, project => $self->{project})) {
      CSG::Storage::Slots::Exceptions::SlotExists->throw();
    }

    my $pool = $project->next_available_pool(
      parent  => $self->parent,
      size    => $self->{size},
    );

    unless ($pool) {
      CSG::Storage::Slots::Exceptions::Pools::NoPoolAvailable->throw();
    }

    my $record = $schema->resultset('Slot')->find_or_create(
      {
        name    => $self->{name},
        size    => $self->{size},
        sha1    => sha1_hex($self->{name}),
        pool_id => $pool->id,
      }
    );

    $self->_record($record);
  }

  return $self->$orig(@_);
};

sub _set_name {
  my ($self, $new, $old) = @_;

  if ($new =~ /^([^-]+)\-.*$/) {
    $self->parent($1);

    my $slot = $self->find(
      prefix  => $self->prefix,
      name    => $self->parent,
      project => $self->project,
      prefix  => $self->prefix,
    );

    unless ($slot) {
      CSG::Storage::Slots::Exceptions::Slot::Parent::DoesNotExist->throw();
    }
  }
}

sub _set_size {
  my ($self, $new, $old) = @_;

  if ($self->has_record) {
    $self->_record->update({size => $new});
  }
}

sub _set_project {
  my ($self, $new, $old) = @_;

  my $schema = CSG::Storage::Slots::DB->new();

  unless ($schema->resultset('Project')->find({name => $new})) {
    CSG::Storage::Slots::Exceptions::Project::DoesNotExist->throw();
  }
}

sub _build_path {
  my ($self) = @_;
  my $pool = $self->_record->pool;
  return File::Spec->join($self->prefix, $pool->hostname, $pool->path, (split(//, $self->sha1))[0 .. 3], $self->name);
}

sub to_string {
  return shift->path;
}

sub is_empty {
  return (Path::Tiny::path(shift->path)->children) ? $FALSE : $TRUE;
}

sub find {
  my $class  = shift;
  my %params = @_;
  my $schema = CSG::Storage::Slots::DB->new();
  my $slot   = $schema->resultset('Slot')->find_slot($params{name}, $params{project});

  return unless $slot;

  return $class->new(
    prefix  => $params{prefix},
    name    => $slot->name,
    project => $slot->pool->project->name,
    size    => $slot->size,
    _record => $slot,
  );
}

sub find_or_create {
  my $class = shift;
  return $class->find(@_) // $class->new(@_);
}

sub exists {
  my $class  = shift;
  my %params = @_;
  my $schema = CSG::Storage::Slots::DB->new();
  return defined $schema->resultset('Slot')->find_slot($params{name}, $params{project});
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
