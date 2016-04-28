package CSG::Storage::SlotFS::Command::init;

use CSG::Storage::SlotFS -command;

use Modern::Perl;
use Module::Load;
use Try::Tiny;
use Number::Bytes::Human qw(parse_bytes);

use CSG::Logger;

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

    my %params = (
      name => $opts->{name},
      size => parse_bytes($opts->{size}),
    );

    if ($opts->{prefix}) {
      $params{prefix} = $opts->{prefix};
    }

    $slot = $class->new(%params);

    $slot->initialize;
    $slot->save_manifest;

    $logger->info($slot->to_string);
  }
  catch {
    if (not ref $_) {
      $logger->error('Unknown error occured');
    } elsif ($_->isa('CSG::Storage::Slots::Exceptions::Sample::FailedSkeletonDirectory')) {
      $logger->error($_->error);
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

CSG::Storage::SlotFS::Command::init - Initialize slot directory structure
