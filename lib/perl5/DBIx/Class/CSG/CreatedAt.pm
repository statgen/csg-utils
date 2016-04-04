package DBIx::Class::CSG::CreatedAt;

use base qw(DBIx::Class);

use Modern::Perl;
use DateTime;
use DateTime::Format::MySQL;

sub insert {
  my $self = shift;
  my $rec  = $self->next::method(@_);
  my $now  = DateTime->now(time_zone => 'America/Detroit');

  $rec->update(
    {
      created_at => DateTime::Format::MySQL->format_datetime($now),
    }
  );

  return $rec;
}

1;
