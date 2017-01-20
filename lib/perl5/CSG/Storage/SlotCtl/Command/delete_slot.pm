package CSG::Storage::SlotCtl::Command::delete_slot;

use CSG::Storage::SlotCtl -command;
use CSG::Base qw(file);
use CSG::Constants;
use CSG::Storage::Slots::DB;

use IO::Prompter {ask => [-in => *STDIN, -out => *STDOUT, -yn, -style => 'bold red']};

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['name|n=s', 'Name of slot to delete', {required => 1}],
    ['force|f',  'Delete slot and contents'],
    ['yes|y',    'Answer yes to prompts'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $project = $schema->resultset('Project')->find({name => $self->app->global_options->{project}});

  unless ($project) {
    $self->exit_with_error('project does not exist');
  }

  unless ($project->has_slot($opts->{name})) {
    $self->exit_with_error('slot does not exist');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $slot = $schema->resultset('Slot')->find_slot($opts->{name}, $self->app->global_options->{project});
  my $slotfs = CSG::Storage::Slots->find(
    name    => $slot->name,
    project => $slot->pool->project->name,
    prefix  => $self->app->global_options->{prefix},
  );

  if ($self->app->global_options->{verbose}) {
    say "Slot Details:\t" . $slot->to_string;
  }

  unless ($opts->{yes}) {
    exit 1 unless ask "Really delete slot, $opts->{name}, and all it's contents? [yn]";
  }

  if (not $slotfs->is_empty and not $opts->{force}) {
    say 'Slot is not empty, refusing to delete';
    exit 1;
  }

  my $deleted = remove_tree($slotfs->path);
  $slot->delete();

  if ($self->app->global_options->{verbose}) {
    say "Deleted slot record and removed $deleted files from " . $slotfs->path;
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::delete_slot - Delete an existing slot. Use with caution.
