package CSG::Storage::SlotCtl::Command::add_type;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (['name|n=s', 'Name for type', {required => 1}],);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($schema->resultset('Type')->find({name => $opts->{name}})) {
    $self->usage_error('type already exists');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;
  $schema->resultset('Type')->create({name => $opts->{name}});
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::add_type - Add a new pool type
