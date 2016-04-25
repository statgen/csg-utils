package CSG::Storage::SlotCtl::Command::update_slot;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots::DB;
use CSG::Base;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['id=i', 'Slot id to modifiy', {required => 1}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;
}

sub execute {
  my ($self, $opts, $args) = @_;

  # TODO - allow the size to be changed but nothing else
  #        if the size does change then need to apply the
  #        slot allocation logic to make sure we don't
  #        go over the available space. if we do bail out
  #        without making changes.
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::update_slot - Update an existing slot
