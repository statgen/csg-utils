package CSG::Mapper::Command::verify;

use CSG::Mapper -command;
use CSG::Base;
use CSG::Constants qw(:mapping);
use CSG::Mapper::DB;

sub opt_spec {
  return (
    ['complete',   'show completed sample results' ],
    ['incomplete', 'show incomplete sample results'],
  );
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $build   = $self->app->global_options->{build};
  my $schema  = CSG::Mapper::DB->new();
  my $results = $schema->resultset('Result')->search(
    {
      'me.build'   => $build,
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
      build   => $build,
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
