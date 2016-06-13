package CSG::Storage::SlotCtl::Command::list_projects;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub execute {
  my ($self, $opts, $args) = @_;
  for my $project ($schema->resultset('Project')->all()) {
    say 'Id: ' . $project->id . ' Name: ' . $project->name;
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::list_projects - List all defined projects
