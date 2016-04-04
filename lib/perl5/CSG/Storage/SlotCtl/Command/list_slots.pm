package CSG::Storage::SlotCtl::Command::list_slots;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots::DB;

use Modern::Perl;

sub opt_spec {
  return (
    ['project|p=s', 'Limit slots displayed to a specific project'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{project}) {
    my $schema = CSG::Storage::Slots::DB->new();
    unless ($schema->resultset('Project')->find({name => $opts->{project}})) {
      $self->usage_error('Project does not exist');
    }
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Storage::Slots::DB->new();
  my $slots;

  if ($opts->{project}) {
    my $project = $schema->resultset('Project')->find({name => $opts->{project}});
    for my $pool ($project->pools) {
      for my $slot ($pool->slots) {
        say $slot->to_string;
      }
    }
  } else {
    for my $slot ($schema->resultset('Slot')->all()) {
      say $slot->to_string;
    }
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::list_slots - List slots for a project or all slots
