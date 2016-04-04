## no critic (RequireArgUnpacking, ProhibitNestedSubs, RequireFilenameMatchesPackage, ProhibitMultiplePackages)
package CSG::Mapper::Logger;

use CSG::Base qw(logging);
use CSG::Constants;
use CSG::Mapper::Config;
use CSG::Mapper::Logger::Dispatch::DBI;

use Moose;

has 'job_id' => (
  is        => 'rw',
  isa       => 'Int',
  predicate => 'has_job_id',
  trigger   => \&_set_job_id,
);

has '_logger' => (
  is      => 'rw',
  isa     => 'Log::Dispatch',
  lazy    => 1,
  builder => '_build_logger',
);

sub _set_job_id {
  my ($self, $new, $prev) = @_;

  if ($self->has_job_id) {
    my $conf = CSG::Mapper::Config->new(project => 'topmed');

    $self->_logger->add(
      CSG::Mapper::Logger::Dispatch::DBI->new(
        dbh       => DBI->connect($conf->dsn, $conf->get('db', 'user'), $conf->get('db', 'pass'), {mysql_auto_reconnect => 1}),
        table     => 'logs',
        min_level => 'info',
      )
    );
  }

  $self->{job_id} = $new;

  return;
}

sub _build_logger {
  my ($self) = @_;

  sub _add_timestamp {
    my (%log) = @_;

    my $timestamp = DateTime->now(time_zone => $TIMEZONE);
    my $level = uc($log{level});

    return ($log{job_id})
      ? qq($timestamp [$level] job[$log{job_id}] $log{message})
      : qq($timestamp [$level] $log{message});
  }

  my $log = Log::Dispatch->new();

  $log->add(
    Log::Dispatch::Screen->new(
      stdout    => 1,
      stderr    => 0,
      newline   => 1,
      min_level => 'debug',
      max_level => 'warning',
      callbacks => \&_add_timestamp,
    )
  );

  $log->add(
    Log::Dispatch::Screen->new(
      stdout    => 0,
      stderr    => 1,
      newline   => 1,
      min_level => 'error',
      max_level => 'emergency',
      callbacks => \&_add_timestamp,
    )
  );

  return $log;
}

sub _log {
  my ($self, $level, $msg) = @_;
  return $self->_logger->log(level => $level, message => $msg, job_id => $self->job_id);
}

sub debug     {return shift->_log('debug',     @_);}
sub info      {return shift->_log('info',      @_);}
sub notice    {return shift->_log('notice',    @_);}
sub warning   {return shift->_log('warning',   @_);}
sub error     {return shift->_log('error',     @_);}
sub critical  {return shift->_log('critical',  @_);}
sub alert     {return shift->_log('alert',     @_);}
sub emergency {return shift->_log('emergency', @_);}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
