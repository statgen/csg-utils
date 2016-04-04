package CSG::Mapper::Command::export;

use CSG::Mapper -command;
use CSG::Base qw(cmd);
use CSG::Mapper::Config;
use CSG::Mapper::DB;

sub opt_spec {
  return (['build=s', 'Reference build used to mapping these samples'],);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Mapper::DB->new();
  my $logger = CSG::Mapper::Logger->new();

  $self->{stash}->{schema} = $schema;
  $self->{stash}->{logger} = $logger;

  unless ($self->app->global_options->{project}) {
    $self->usage_error('Project is required');
  }

  my $config = CSG::Mapper::Config->new(project => $self->app->global_otions->{project});
  $self->{stash}->{config} = $config;

  unless ($self->can('_export_' . $self->app->global_options->{project})) {
    $self->usage_error('No export method defined for this project');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $project = $self->app->global_options->{project};
  my $schema  = $self->{stash}->{schema};
  my $results = $schema->resultset('Result')->search(
    {
      'me.build'           => $opts->{build},
      'sample.exported_at' => undef,
      'state.name'         => 'completed',
      'project.name'       => $project,
    }, {
      join => ['state', {sample => 'project'}],
    },
  );

  while (my $result = $results->next) {
    my $export_meth = qq{_export_$project};
    $self->$export_meth($result->sample, $result->build);
  }
}

sub _export_topmed {
  my ($self, $sample, $build) = @_;

  my $logger = $self->{stash}->{logger};
  my $cmd = sprintf '/usr/cluster/monitor/bin/topmedcmd.pl %s mapped%d completed', $sample->sample_id, $build;

  $logger->debug("EXPORT CMD: '$cmd'") if $self->app->global_options->{debug};

  try {
    run($cmd) unless $self->app->global_options->{dry_run};
    $logger->info('Exported sample[' . $sample->sample_id . '] to topmedcmd ') if $self->app->global_options->{verbose};
    $sample->update({exported_at => DateTime->now()});
  }
  catch {
    unless (ref $_) {
      $logger->error($_);
    }
  };
}

1;

__END__

=head1

CSG::Mapper::Command::export - export remapping jobs
