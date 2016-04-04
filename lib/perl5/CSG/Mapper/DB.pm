package CSG::Mapper::DB;

use base qw(CSG::Mapper::DB::Schema);

use CSG::Base;
use CSG::Constants;
use CSG::Mapper::Config;
use CSG::Mapper::DB::Schema;

sub new {
  my $conf = CSG::Mapper::Config->new(project => 'topmed');
  return __PACKAGE__->connect($conf->dsn, $conf->get('db', 'user'), $conf->get('db', 'pass'));
}

sub now {
  my $now = DateTime->now(time_zone => $TIMEZONE);
  return DateTime::Format::MySQL->format_datetime($now),;
}

1;
