package CSG::Storage::SlotFS::Command::path;

use CSG::Storage::SlotFS -command;

use Modern::Perl;
use Module::Load;
use Try::Tiny;

use CSG::Logger;
use CSG::Storage::Slots::DB;

sub opt_spec {
  return (
    ['project|p=s', 'Project the slot belongs to',                              {required => 1}],
    ['name|n=s',    'Slot name to inspect',                                     {required => 1}],
    ['path|t=s',    'Path within the slot to retreive (i.e. incoming, backup)', {required => 1}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Storage::Slots::DB->new();
  my $project = $schema->resultset('Project')->find({name => $opts->{project}});
  unless ($project) {
    $self->usage_error('Project does not exist');
  }

  unless ($schema->resultset('Slot')->find_slot($opts->{name}, $opts->{project})) {
    $self->usage_error('Slot does not exist');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();
  my $class  = 'CSG::Storage::SlotFS::' . ucfirst(lc($opts->{project}));
  my $slot   = undef;
  my $path   = qq{$opts->{path}_path};

  try {
    load $class;
    $slot = $class->new(name => $opts->{name});
    $logger->info($slot->$path);
  }
  catch {
    say $_;
    $rc = 1;
  };

  exit $rc;
}

1;

__END__

=head1

CSG::Storage::SlotFS::Command::path - find various path components for a storage slot
