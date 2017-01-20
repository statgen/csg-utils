package CSG::Storage::SlotCtl::Command::show;

use CSG::Storage::SlotCtl -command;

use CSG::Base qw(formats);
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (['name|n=s', 'Name of slot to display', {required => 1}]);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $name    = $opts->{name};
  my $project = $self->app->global_options->{project};

  unless (CSG::Storage::Slots->exists(name => $name, project => $project)) {
    $self->exit_with_error("Slot, $name, does not exists in project, $project");
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $name    = $opts->{name};
  my $project = $self->app->global_options->{project};
  my $slot    = $schema->resultset('Slot')->find_slot($name, $project);

  say YAML::Dump($slot->to_hashref);
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::show - Display details of a slot
