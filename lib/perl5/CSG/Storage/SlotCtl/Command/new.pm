package CSG::Storage::SlotCtl::Command::new;

use CSG::Storage::SlotCtl -command;
use CSG::Base qw(file);
use CSG::Logger;
use CSG::Storage::Slots;

use Number::Bytes::Human qw(parse_bytes);

sub opt_spec {
  return (
    ['name|n=s', 'Name for the slot',                                        {required => 1}],
    ['size|s=s', 'Disk space required in human readable format (i.e. 400G)', {required => 1}],
    # TODO - need option to disable pool host exclusion
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  my $name    = $opts->{name};
  my $project = $self->app->global_options->{project};

  if (CSG::Storage::Slots->exists(name => $name, project => $project)) {
    $self->usage_error("Slot, $name, already exists in project, $project");
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();
  my $slot   = undef;

  try {
    $slot = CSG::Storage::Slots->new(
      name    => $opts->{name},
      project => $self->app->global_options->{project},
      size    => parse_bytes($opts->{size}),
      prefix  => $self->app->global_options->{prefix},
    );

    make_path($slot->path) unless -e $slot->path;

    $logger->info($slot->to_string);
  }
  catch {
    if (not ref $_) {
      $logger->error($_);

    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Slot::DoesNotExist')) {
      $logger->critical($_->error);

    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Slot::Parent::DoesNotExist')) {
      $logger->critical($_->error);

    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Project::DoesNotExist')) {
      $logger->critical($_->error);

    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Pool::DoesNotExist')) {
      $logger->critical($_->error);

    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Pools::NoPoolAvailable')) {
      $logger->critical($_->error);

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

CSG::Storage::SlotCtl::Command::new - Create a new storage slot
