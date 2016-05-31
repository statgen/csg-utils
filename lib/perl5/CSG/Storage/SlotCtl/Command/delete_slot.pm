package CSG::Storage::SlotCtl::Command::delete_slot;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

use IO::Prompter {ask => [-in => *STDIN, -out => *STDOUT, -yn, -style => 'bold red']};

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  ## no tidy
  return (
    ['project|p=s', 'Project slot is part of', {required => 1}],
    ['name|n=s',    'Name of slot to delete',  {required => 1}],
  );
  ## use tidy
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $project = $schema->resultset('Project')->find({name => $opts->{project}});
  unless ($project) {
    $self->usage_error('project does not exist');
  }

  unless ($project->has_slot($opts->{name})) {
    $self->usage_error('slot does not exist');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $slot = $schema->resultset('Slot')->find_slot($opts->{name}, $opts->{project});

  exit unless ask "Really delete slot, $opts->{name}? [yn]";

  say 'Slot Details:';
  say $slot->to_string;
  $slot->delete();
  say 'Deleted';
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::delete_slot - Delete an existing slot. Use with caution.
