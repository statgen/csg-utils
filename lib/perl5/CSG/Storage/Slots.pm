package CSG::Storage::Slots;

use CSG::Base;
use CSG::Storage::Slots::DB;
use CSG::Storage::Slots::Exceptions;
use CSG::Types;

use Moose;
use Digest::SHA qw(sha1_hex);
use overload '""' => sub {shift->to_string};

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

has 'exclude' => (
  is      => 'ro',
  isa     => 'Maybe[Int]',
  default => sub {
    return undef;
  },
);

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

has 'prefix' => (
  is      => 'ro',
  isa     => 'ValidPrefixPath',
  default => sub {
    return '/net';
  },
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

    my $pool = $project->next_available_pool(
      size    => $self->{size},
      exclude => $self->exclude,
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

  # TODO - need to keep all associated slots off the same hosts
  #
  if ($new =~ /^([^-]+)\-.*$/) {
    my $name = $1;
    my $slot = $self->find(
      name    => $name,
      project => $self->project,
      prefix  => $self->prefix,
    );

    unless ($slot) {
      CSG::Storage::Slots::Exceptions::Slot::DoesNotExist->throw();
    }

    $self->exclude($slot->pool_id);
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

sub find {
  my $class  = shift;
  my %params = @_;
  my $schema = CSG::Storage::Slots::DB->new();
  my $slot   = $schema->resultset('Slot')->find_slot($params{name}, $params{project});

  return unless $slot;

  return $class->new(
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

no Moose;
__PACKAGE__->meta->make_immutable;

1;
