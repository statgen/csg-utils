package CSG::Mapper;

use App::Cmd::Setup -app;

sub global_opt_spec {
  return (
    ['debug|d',   'Debug output'],
    ['verbose|v', 'Verbose output'],
    ['dry-run|n', 'Dry run; show what would be done without actaully doing anything'],
    ['help|h',    'Usage'],
    ['cluster=s', 'Cluster environment (valid clusters: csg|flux)'],
    ['project=s', 'Project settings to load'],
  );
}

1;
