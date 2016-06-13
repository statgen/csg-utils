package CSG::Storage::SlotCtl::Command::update_project;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (['name=s', 'New project name', {requried => 1}]);
}

sub execute {
  my ($self, $opts, $args) = @_;
  my $project = $schema->resultset('Project')->find({name => $self->app->global_options->{project}});
  $project->update({name => $opts->{name}});
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::update_project - Update a projects details
