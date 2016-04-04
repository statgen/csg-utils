package CSG::Storage::Slots::DB;

use base qw(CSG::Storage::Slots::DB::Schema);

use CSG::Storage::Config;

sub new {
  my $conf = CSG::Storage::Config->new();
  return __PACKAGE__->connect($conf->dsn, $conf->db_user, $conf->db_pass);
}

1;
