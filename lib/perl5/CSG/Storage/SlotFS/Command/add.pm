package CSG::Storage::SlotFS::Command::add;

use CSG::Storage::SlotFS -command;

sub opt_spec {
  return (
    ['project|p=s', 'Project name the slot belongs to',                                         {required => 1}],
    ['name|n=s',    'Name of the slot to initialize',                                           {required => 1}],
    ['size|s=s',    'Size of the initial sample directory in human readable form (i.e. 400GB)', {required => 1}],
    ['prefix|=s',   'Optional path prefix (e.g. /tmp)'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{prefix} and not -e $opts->{prefix}) {
    $self->usage_error('Prefix must exist');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();
  my $class  = 'CSG::Storage::SlotFS::' . ucfirst(lc($opts->{project}));
  my $slot   = undef;

  try {
    load $class;

    (my $name = $opts->{name}) =~ s/^([^-]+)\-.*$/$1/g;

    my %params = (
      name => $name,
      size => parse_bytes($opts->{size}),
    );

    $slot = $class->find(%params);
    unless ($slot) {
      CSG::Storage::Slots::Exceptions::Slot::DoesNotExist->throw();
    }

    $params{exclude} = $slot->pool_id;
    $class->new(%params);
  }
  catch {
    if (not ref $_) {
      $logger->error('Unknown error occured');
    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Sample::FailedSkeletonDirectory')) {
      $logger->error($_->description);
    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Slot::DoesNotExist')) {
      $logger->error($_->description);
    } else {
      if ($_->isa('Exception::Class')) {
        $logger->error($_->error);
      } else {
        $logger->error("something went sideways: $_");
      }
    }

    $rc = 1;
  }
  finally {
    unless (@_) {
      $logger->info($slot->to_string);
    }
  };

  exit $rc;
}

1;

__END__

=head1

CSG::Storage::SlotFS::Command::add - Add a slot to an existing slot directory
