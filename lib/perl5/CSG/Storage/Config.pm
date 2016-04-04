package CSG::Storage::Config;

use FindBin qw($Bin);
use Modern::Perl;
use Config::Tiny;
use Moose;

use CSG::Storage::Types;

has '_file'    => (is => 'ro', isa => 'ValidFile',    lazy => 1, builder => '_build_conf_file');
has '_content' => (is => 'ro', isa => 'Config::Tiny', lazy => 1, builder => '_build_content');

has 'db'      => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_db');
has 'db_host' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_db_host');
has 'db_user' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_db_user');
has 'db_pass' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_db_pass');

sub _build_conf_file {
  return $ENV{SLOTS_CONF} // qq{$Bin/../config/slots.ini};
}

sub _build_content {
  return Config::Tiny->read(shift->_file);
}

sub _build_db {
  return $ENV{SLOTS_DB} // shift->get('database', 'db');
}

sub _build_db_host {
  return $ENV{SLOTS_DB_HOST} // shift->get('database', 'host');
}

sub _build_db_user {
  return $ENV{SLOTS_DB_USER} // shift->get('database', 'username');
}

sub _build_db_pass {
  return $ENV{SLOTS_DB_PASS} // shift->get('database', 'password');
}

sub dsn {
  my ($self) = @_;
  return sprintf 'dbi:mysql:database=%s;host=%s;port=3306', $self->db, $self->db_host;
}

sub get {
  my ($self, $category, $name) = @_;
  return $self->_content->{$category}->{$name};
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
