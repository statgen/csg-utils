package CSG::Mapper::Command::log;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Mapper::Logger;

sub opt_spec {
  return (
    ['meta-id=i', 'Job meta id (database id)'],
    ['message=s', 'Text of message to log', {required => 1}],
    [
      'level=s',
      'Log level for this message (valid levels: debug|info|notice|warning|error|critical|alert|emergency) [default: info]', {
        default   => 'info',
        callbacks => {
          regex => sub {
            shift =~ /debug|info|notice|warning|error|critical|alert|emergency/;
          }
        }
      }
    ],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{meta_id}) {
    my $schema = CSG::Mapper::DB->new();
    unless ($schema->resultset('Job')->find($opts->{meta_id})) {
      $self->usage_error("Job id, $opts->{meta_id}, does not exist!");
    }
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $logger = CSG::Mapper::Logger->new();
  if ($opts->{meta_id}) {
    $logger->job_id($opts->{meta_id});
  }

  my $level = $opts->{level};
  $logger->$level($opts->{message});
}

1;

__END__

=head1

CSG::Mapper::Command::log - log remapping job info
