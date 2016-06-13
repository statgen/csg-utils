package CSG::Storage::SlotCtl::Command::list_pools;

use CSG::Storage::SlotCtl -command;
use CSG::Base;
use CSG::Storage::Slots::DB;

my $schema = CSG::Storage::Slots::DB->new();

sub execute {
  my ($self, $opts, $args) = @_;

  my $pools = $schema->resultset('Pool');
  if ($opts->{project}) {
    $pools = $pools->search({'project.name' => $self->app->global_options->{project}}, {join => 'project'});
  }

  for my $pool ($pools->all()) {
    say 'Name: ' . $pool->name;
    say "\tID: " . $pool->id;
    say "\tProject: " . $pool->project->name;
    say "\tHostname: " . $pool->hostname;
    say "\tPath: " . $pool->path;
    say "\tSpace Total: " . $pool->size_total;
    say "\tSpace Used: " . $pool->size_used;
    say '';
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::list_pools - List all defined pools
