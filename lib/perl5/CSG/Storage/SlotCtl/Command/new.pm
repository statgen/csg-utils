package CSG::Storage::SlotCtl::Command::new;

use CSG::Storage::SlotCtl -command;
use CSG::Storage::Slots;
use CSG::Logger;

use Modern::Perl;
use Try::Tiny;
use Number::Bytes::Human qw(parse_bytes);

sub opt_spec {
  return (
    ['name|n=s',    'Name for the slot',                                        {required => 1}],
    ['project|p=s', 'Project name for slot filesystems [default: topmed]',      {default  => 'topmed'}],
    ['size|s=s',    'Disk space required in human readable format (i.e. 400G)', {required => 1}],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $rc     = 0;
  my $logger = CSG::Logger->new();
  my $slot   = undef;

  try {
    $slot = CSG::Storage::Slots->new(
      name    => $opts->{name},
      project => $opts->{project},
      size    => parse_bytes($opts->{size}),
    );
  }
  catch {
    if (not ref $_) {
      $logger->error($_);
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

CSG::Storage::SlotCtl::Command::new - Create a new storage slot
