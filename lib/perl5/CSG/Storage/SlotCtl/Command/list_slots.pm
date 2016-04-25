package CSG::Storage::SlotCtl::Command::list_slots;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots::DB;
use CSG::Base;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (['project|p=s', 'Limit slots displayed to a specific project'],);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{project}) {
    unless ($schema->resultset('Project')->find({name => $opts->{project}})) {
      $self->usage_error('Project does not exist');
    }
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $slots = $schema->resultset('Slot');
  if ($opts->{project}) {
    $slots = $slots->search(
      {
        'project.name' => $opts->{project},
      },
      {
        join => {pool => 'project'}
      }
    );
  }

  for my $slot ($slots->all()) {
    say $slot->to_string;
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::list_slots - List slots for a project or all slots
