package CSG::Mapper;

use App::Cmd::Setup -app;
use CSG::Base;
use CSG::Constants qw(:basic :mapping);

sub global_opt_spec {
  return (
    ['help|h',    'Usage'],
    ['debug|d',   'Debug output'],
    ['verbose|v', 'Verbose output'],
    ['dry-run|n', 'Dry run; show what would be done without actaully doing anything'],
    [
      'cluster=s',
      'Cluster environment (valid clusters: csg|flux) [default: csg]', {
        default   => 'csg',
        callbacks => {
          regex => sub {shift =~ /$VALID_CLUSTER_REGEXPS/}
        }
      }
    ],
    ['project=s', 'Project settings to load [default: topmed]', {default => 'topmed'}],
    [
      'build=i',
      'Refernce build (valid values: 37 or 38) [default: 38]', {
        default   => '38',
        callbacks => {
          regex => sub {
            shift =~ /37|38/
          }
        }
      }
    ],
  );
}

1;
