package CSG::Mapper::Logger::Dispatch::DBI;

use base qw(Log::Dispatch::DBI);
use CSG::Base;

sub create_statement {
  my $self = shift;
  my $sql  = qq{insert into $self->{table} (job_id, level, message) values (?, ?, ?)};
  return $self->{dbh}->prepare($sql);
}

sub log_message {
  my $self   = shift;
  my %params = @_;
  return $self->{sth}->execute(@params{(qw(job_id level message))});
}

1;
