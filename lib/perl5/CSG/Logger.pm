package CSG::Logger;

use Log::Dispatch;
use Log::Dispatch::Screen;

use Moose;

has '_logger' => (
  is      => 'rw',
  isa     => 'Log::Dispatch',
  lazy    => 1,
  builder => '_build_logger',
  handles => [
    qw(
      debug
      info
      notice
      warning
      error
      critical
      alert
      emergency
    )
  ],
);

sub _build_logger {
  my ($self) = @_;

  my $log = Log::Dispatch->new();

  $log->add(
    Log::Dispatch::Screen->new(
      stdout    => 1,
      stderr    => 0,
      newline   => 1,
      min_level => 'debug',
      max_level => 'warning',
    )
  );

  $log->add(
    Log::Dispatch::Screen->new(
      stdout    => 0,
      stderr    => 1,
      newline   => 1,
      min_level => 'error',
      max_level => 'emergency',
    )
  );

  return $log;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
