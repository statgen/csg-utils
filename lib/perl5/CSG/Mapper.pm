package CSG::Mapper;

use App::Cmd::Setup -app;

use CSG::Base;
use CSG::Mapper::Util qw(detect_cluster);

sub global_opt_spec {
  return (
    ['help|h',    'Usage'],
    ['debug|d',   'Debug output'],
    ['verbose|v', 'Verbose output'],
    ['dry-run|n', 'Dry run; show what would be done without actaully doing anything'],
    ['cluster=s', 'Cluster environment (valid clusters: csg|flux)'                  ],
    ['project=s', 'Project settings to load'                                        ],
  );
}

sub prepare_args {
  my ($self) = @_;

  if (my $cluster = detect_cluster()) {
    unless (grep {/\-\-cluster/} @ARGV) {
      push @ARGV, ('--cluster', $cluster);
    }
  }

  return (@ARGV);
}

sub validate_args {
  my ($self, $opt, $args) = @_;
  die;
}

1;
