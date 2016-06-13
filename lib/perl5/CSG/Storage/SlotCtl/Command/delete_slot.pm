package CSG::Storage::SlotCtl::Command::delete_slot;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

use IO::Prompter {ask => [-in => *STDIN, -out => *STDOUT, -yn, -style => 'bold red']};

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (['name|n=s', 'Name of slot to delete', {required => 1}]);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $project = $schema->resultset('Project')->find({name => $self->app->global_options->{project}});

  unless ($schema->resultset('Project')->find({name => $self->app->global_options->{project}})) {
    $self->usage_error('project does not exist');
  }

  unless ($project->has_slot($opts->{name})) {
    $self->usage_error('slot does not exist');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $slot = $schema->resultset('Slot')->find_slot($opts->{name}, $self->app->global_options->{project});

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
