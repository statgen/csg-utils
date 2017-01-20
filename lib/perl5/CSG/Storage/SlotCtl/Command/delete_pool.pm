package CSG::Storage::SlotCtl::Command::delete_pool;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

use IO::Prompter {ask => [-in => *STDIN, -out => *STDOUT, -yn, -style => 'bold red']};

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['pool-id|i=i', 'ID of the pool to delete', {required => 1}]
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($schema->resultset('Pool')->find($opts->{pool_id})) {
    $self->exit_with_error("Pool id, $opts->{pool_id}, does not exist!");
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $pool = $schema->resultset('Pool')->find($opts->{pool_id});

  die 'Slots are still allocated to this pool!' if $pool->slots->count > 0;

  say 'Pool Details:';
  say 'Name: '              . $pool->name;
  say "\tID: "              . $pool->id;
  say "\tProject: "         . $pool->project->name;
  say "\tHostname: "        . $pool->hostname;
  say "\tPath: "            . $pool->path;
  say "\tSpace Total: "     . $pool->size_total;
  say "\tSpace Used: "      . $pool->size_used;
  say "\tSlots allocated: " . $pool->slots->count;
  say '';

  exit unless ask 'Really delete pool, ' . $pool->name . '? [yn]';
  $pool->delete();
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::delete_pool - Delete a pool. Use with caution.
