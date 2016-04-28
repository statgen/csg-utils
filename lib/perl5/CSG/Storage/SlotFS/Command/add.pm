package CSG::Storage::SlotFS::Command::add;

use CSG::Storage::SlotFS -command;
use CSG::Logger;
use CSG::Base;

use Module::Load;
use Number::Bytes::Human qw(parse_bytes);

sub opt_spec {
  return (
    ['project|p=s', 'Project name the slot belongs to',                                         {required => 1}],
    ['name|n=s',    'Name of the slot to initialize',                                           {required => 1}],
    ['size|s=s',    'Size of the initial sample directory in human readable form (i.e. 400GB)', {required => 1}],
    ['prefix=s',    'Optional path prefix (e.g. /tmp)',                                         {default  => '/net'}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless (-e $opts->{prefix}) {
    $self->usage_error('Prefix must exist');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();

  try {
    my $class = 'CSG::Storage::SlotFS::' . ucfirst(lc($opts->{project}));

    load $class;

    (my $name = $opts->{name}) =~ s/^([^-]+)\-.*$/$1/g;

    my $slot = $class->find(
      name    => $name,
      project => $opts->{project},
      prefix  => $opts->{prefix},
    );

    unless ($slot) {
      CSG::Storage::Slots::Exceptions::Slot::DoesNotExist->throw();
    }

    my $new_slot = $class->new(
      name    => $opts->{name},
      size    => parse_bytes($opts->{size}),
      exclude => $slot->pool_id,
      project => $opts->{project},
      prefix  => $opts->{prefix},
    );

    $logger->info($new_slot->to_string);
    $new_slot->save_manifest;
  }
  catch {
    if (not ref $_) {
      $logger->error($_);
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
  };

  exit $rc;
}

1;

__END__

=head1

CSG::Storage::SlotFS::Command::add - Add a slot to an existing slot directory
