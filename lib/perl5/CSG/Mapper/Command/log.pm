package CSG::Mapper::Command::log;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Mapper::DB;
use CSG::Mapper::Logger;

sub opt_spec {
  return (
    ['meta-id=i', 'Job meta id (database id)'],
    ['message=s', 'Text of message to log'],
    ['level=s',   'Log level for this message (valid levels:debug|info|notice|warning|error|critical|alert|emergency)'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $logger = CSG::Mapper::Logger->new();
  my $level = $opts->{level} // q{info};

  unless ($logger->can($level)) {
    $self->usage_error('Invalid loglevel');
  }

  if ($opts->{meta_id}) {
    $logger->job_id($opts->{meta_id});
  }

  $self->{stash}->{logger} = $logger;
  $self->{stash}->{level}  = $level;
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $logger = $self->{stash}->{logger};
  my $level  = $self->{stash}->{level};

  $logger->$level($opts->{message});
}

1;

__END__

=head1

CSG::Mapper::Command::log - log remapping job info
