package CSG::Storage::SlotCtl::Command::add_project;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots::DB;

use Modern::Perl;

sub opt_spec {
  return (['name|n=s', 'Project name', {required => 1}]);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Storage::Slots::DB->new();
  if ($schema->resultset('Project')->find({name => $opts->{project}})) {
    $self->usage_error('project already exists');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Storage::Slots::DB->new();
  $schema->resultset('Project')->create({name => $opts->{name}});
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::add_project - Add a new project
