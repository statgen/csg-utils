package App::Cmd::Plugin::Logger;

{
  $App::Cmd::Plugin::Logger::VERSION = '0.1';
}

use App::Cmd::Setup -plugin => {
  exports => [qw(exit_with_error)],
};

use CSG::Base;
use CSG::Logger;

sub exit_with_error {
  my ($plugin, $cmd, $app, $msg) = @_;
  my $logger = CSG::Logger->new();
  $logger->error($msg);
  exit 1;
}

1;
