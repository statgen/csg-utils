package CSG::Storage::SlotCtl::Command::list_slots;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub execute {
  my ($self, $opts, $args) = @_;

  my $slots = $schema->resultset('Slot')->search(
    {
      'project.name' => $self->app->global_options->{project},
    }, {
      join => {pool => 'project'}
    }
  );

  for my $slot ($slots->all()) {
    say $slot->to_string;
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::list_slots - List slots for a project or all slots
