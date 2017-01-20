package CSG::Storage::SlotCtl::Command::update_slot;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

use Number::Bytes::Human qw(parse_bytes);

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['id=i',   'Slot id to modifiy',                                          {required => 1}],
    ['size=s', 'New size for the slot in human readable format (e.g. 200GB)', {required => 1}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($schema->resultset('Slot')->find($opts->{id})) {
    $self->exit_with_error('Unable to locate slot');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $size = parse_bytes($opts->{size});
  my $slot = $schema->resultset('Slot')->find($opts->{id});

  if ($size > $slot->size) {
    unless ($slot->pool->is_available($size - $slot->size)) {
      say 'Pool does not have enough available space';
      exit 1;
    }
  }

  $slot->update({size => $size});
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::update_slot - Update the size of an existing slot
