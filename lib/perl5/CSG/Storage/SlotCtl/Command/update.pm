package CSG::Storage::SlotCtl::Command::update;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Logger;
use CSG::Storage::Slots;

sub opt_spec {
  return (['name|n=s', 'Slot name', {required => 1}], ['size|s=s', 'New size for slot', {required => 1}]);
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();
  my $slot   = CSG::Storage::Slots->find(
    name    => $opts->{name},
    project => $self->app->global_options->{project}
  );

  if ($slot) {
    $slot->size($opts->{size});
  } else {
    $logger->error('slot not found');
    $rc = 1;
  }

  exit $rc;
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::update - Update an existing storage slot size
