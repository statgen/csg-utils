package CSG::Mapper::Config;

use Moose;

use CSG::Base qw(config file);
use CSG::Constants;
use CSG::Types;

has 'project' => (
  is       => 'ro',
  isa      => 'Str',
  required => 1,
  trigger  => \&_set_project,
);

has '_config_dir' => (
  is      => 'ro',
  isa     => 'Directory',
  default => sub { return dirname($ENV{CSG_MAPPING_CONF}); },
);

has '_global_conf_file' => (
  is      => 'ro',
  isa     => 'FileOnDisk',
  lazy    => 1,
  builder => '_build_global_conf_file',
);

has '_project_conf_file' => (
  is      => 'ro',
  isa     => 'FileOnDisk',
  lazy    => 1,
  builder => '_build_project_conf_file',
);

has '_global_conf' => (
  is      => 'ro',
  isa     => 'Config::Tiny',
  lazy    => 1,
  builder => '_build_global_conf',
);

has '_project_conf' => (
  is      => 'ro',
  isa     => 'Config::Tiny',
  lazy    => 1,
  builder => '_build_project_conf',
);

sub _build_global_conf_file {
  return $ENV{CSG_MAPPING_CONF} // File::Spec->join(shift->_config_dir, 'mapper.ini');
}

sub _build_global_conf {
  return Config::Tiny->read(shift->_global_conf_file);
}

sub _build_project_conf_file {
  my ($self) = @_;
  return File::Spec->join($self->_config_dir, $self->project . '.ini');
}

sub _build_project_conf {
  return Config::Tiny->read(shift->_project_conf_file);
}

sub _build_conf {
  return Config::Tiny->read(shift->_global_conf_file);
}

sub _set_project {
  my ($self, $name, $old_name) = @_;
  croak "no configuraiton for project: $name exists" unless -e $self->_project_conf_file;
}

sub get {
  my ($self, $category, $name) = @_;

  my $section = ($category eq q{global}) ? $UNDERSCORE : $category;

  if (exists $self->_project_conf->{$section} and exists $self->_project_conf->{$section}->{$name}) {
    return $self->_project_conf->{$section}->{$name};
  }

  if (exists $self->_global_conf->{$section} and exists $self->_global_conf->{$section}->{$name}) {
    return $self->_global_conf->{$section}->{$name};
  }

  return;
}

sub dsn {
  return sprintf 'dbi:mysql:database=%s;host=%s;port=%d', @{shift->_global_conf->{db}}{(qw(db host port))};
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
