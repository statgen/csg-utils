package CSG::Storage::SlotCtl::Command::update_pool;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['id=i',    'ID for the slot to modify', {required => 1}],
    ['size=i',  'Total available space in bytes'],
    ['name=s',  'New name for this pool'],
    ['disable', 'Disable the pool'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($schema->resultset('Pool')->find($opts->{id})) {
    $self->usage_error('invalid pool id');
  }

  if ($opts->{name}) {
    if ($schema->resultset('Pool')->find({name => $opts->{name}})) {
      $self->usage_error('pool name already in use');
    }
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $pool = $schema->resultset('Pool')->find($opts->{id});

  if ($opts->{size}) {
    $pool->update({size => $opts->{size}});
  }

  if ($opts->{name}) {
    $pool->update({name => $opts->{name}});
  }

  if ($opts->{disable}) {
    $pool->update({active => 0});
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::update_pool - Update a pools settings
