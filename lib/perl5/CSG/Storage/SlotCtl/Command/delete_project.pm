package CSG::Storage::SlotCtl::Command::delete_project;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

use IO::Prompter {ask => [-in => *STDIN, -out => *STDOUT, -yn, -style => 'bold red']};

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['project-id|i=i', 'ID of project to delete', {required => 1}]
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($schema->resultset('Project')->find($opts->{project_id})) {
    $self->exit_with_error("Project id, $opts->{project_id}, does not exist!");
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $project = $schema->resultset('Project')->find($opts->{project_id});

  die 'Pools are still allocated to this project!' if $project->pools->count > 0;

  say 'Project Details:';
  say 'Id: ' . $project->id . ' Name: ' . $project->name;

  exit unless ask 'Really delete project, ' . $project->name . '? [yn]';
  $project->delete();
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::delete_project - Delete a project. Use with caution.
