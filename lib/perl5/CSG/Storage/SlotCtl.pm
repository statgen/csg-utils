package CSG::Storage::SlotCtl;

use App::Cmd::Setup -app => {plugins => ['Logger']};
use CSG::Base;
use CSG::Constants qw(:basic :mapping);

sub global_opt_spec {
  return (
    ['debug|d',   'Debug output'],
    ['verbose|v', 'Verbose output'],
    ['dry-run',   'Dry run; show what would be done without actaully doing anything'],
    ['project=s', 'Project settings to load [default: topmed]', {default => 'topmed'}],
    [
      'prefix=s',
      'Path prefix for pools [default: /net]', {
        default   => '/net',
        callbacks => {
          regex => sub {shift =~ /\/net|\/dept\/csg|\/tmp/}
        }
      }
    ], [
      'cluster=s',
      'Cluster environment (valid clusters: csg|flux) [default: csg]', {
        default   => 'csg',
        callbacks => {
          regex => sub {
            shift =~ /$VALID_CLUSTER_REGEXPS/;
          }
        }
      }
    ],

  );
}

1;
