package CSG::Storage::SlotCtl::Command::find;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Logger;
use CSG::Storage::Slots;

sub opt_spec {
  return (['name|n=s', 'Slot name', {required => 1}]);
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();
  my $slot   = CSG::Storage::Slots->find(
    name    => $opts->{name},
    project => $self->app->global_options->{project},
    prefix  => $self->app->global_options->{prefix},
  );

  if ($slot) {
    $logger->info($slot->to_string);
  } else {
    $logger->error('slot not found');
    $rc = 1;
  }

  exit $rc;
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::find - Find an existing storage slot path
