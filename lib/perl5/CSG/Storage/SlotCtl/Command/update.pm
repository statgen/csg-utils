package CSG::Storage::SlotCtl::Command::update;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots;
use CSG::Logger;

use Modern::Perl;

sub opt_spec {
  return (
    ['name|n=s',    'Slot name',         {required => 1}],
    ['project|p=s', 'Project name',      {required => 1}],
    ['size|s=s',    'New size for slot', {required => 1}],
  );
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();
  my $slot   = CSG::Storage::Slots->find(name => $opts->{name}, project => $opts->{project});

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
