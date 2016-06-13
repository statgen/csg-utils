package CSG::Storage::SlotCtl::Command::add_project;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (['name|n=s', 'Project name', {required => 1}]);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($schema->resultset('Project')->find({name => $self->app->global_options->{project}})) {
    $self->usage_error('project already exists');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;
  $schema->resultset('Project')->create({name => $opts->{name}});
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::add_project - Add a new project
