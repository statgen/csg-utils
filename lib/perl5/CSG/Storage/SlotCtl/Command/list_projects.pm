package CSG::Storage::SlotCtl::Command::list_projects;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots::DB;

use Modern::Perl;

sub execute {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Storage::Slots::DB->new();
  for my $project ($schema->resultset('Project')->all()) {
    say 'Id: ' . $project->id . ' Name: ' . $project->name;
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::list_projects - List all defined projects
