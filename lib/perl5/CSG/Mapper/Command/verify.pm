package CSG::Mapper::Command::verify;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Constants qw(:mapping);
use CSG::Mapper::DB;

sub opt_spec {
  return (
    ['build=s', 'reference build to verify results'],
    ['complete', 'show completed sample results'],
    ['incomplete', 'show incomplete sample results'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($opts->{build}) {
    $self->usage_error('build is required');
  }

  unless ($opts->{build} =~ /37|38/) {
    $self->usage_error('invalid build');
  }

  unless ($self->app->global_options->{cluster}) {
    $self->usage_error('Cluster environment is required');
  }

  unless ($self->app->global_options->{cluster} =~ /$VALID_CLUSTER_REGEXPS/) {
    $self->usage_error('Invalid cluster environment');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $schema  = CSG::Mapper::DB->new();
  my $results = $schema->resultset('Result')->search(
    {
      'me.build'   => $opts->{build},
      'state.name' => 'completed',
    },
    {
      join => 'state',
    }
  );

  for my $result ($results->all()) {
    my $sample = CSG::Mapper::Sample->new(
      cluster => $self->app->global_options->{cluster},
      record  => $result->sample,
      build   => $opts->{build},
    );

    if ($opts->{complete}) {
      next unless $sample->is_complete;
      say $result->status_line . 'Complete: YES';

    } elsif ($opts->{incomplete}) {
      next if $sample->is_complete;
      say $result->status_line . 'Complete: NO';

    } else {
      say $result->status_line . 'Complete: ' . ($sample->is_complete) ? 'YES' : 'NO';
    }

  }
}

1;

__END__

=head1

CSG::Mapper::Command::verify - verify remapping jobs
