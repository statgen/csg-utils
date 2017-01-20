package CSG::Storage::SlotCtl::Command::collect;

use CSG::Storage::SlotCtl -command;
use CSG::Base qw(file);
use CSG::Logger;
use CSG::Storage::Slots::DB;

use Filesys::DiskUsage qw(du);

my $schema = CSG::Storage::Slots::DB->new();

sub opt_spec {
  return (
    ['prefix=s', 'PREFIX path to apply to NFS filesystems'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{prefix} and not -e $opts->{prefix}) {
    $self->exit_with_error('PREFIX does not exist');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $prefix = $opts->{prefix} // '/net';
  my $logger = CSG::Logger->new();

  for my $fs ($schema->resultset('Pool')->all()) {
    my $path = File::Spec->canonpath(File::Spec->join($prefix, $fs->hostname, $fs->path));

    unless (-e $path) {
      $logger->error("path, $path, does not exist");
      next;
    }

    my $used = du($path);
    $logger->info(sprintf '%s[%s] used: %d total: %d', $fs->name, $path, $used, $fs->size_total);
    $fs->update({size_used => $used});
  }
}

1;

__END__

=head1

CSG::Storage::SlotCtl::Command::collect - Collect disk usage info for a slot storage filesystem
