## no critic (ProhibitPostfixControls, ProhibitNegativeExpressionsInUnlessAndUntilConditions)
package CSG::Mapper::Sample;

use CSG::Base qw(file);
use CSG::Mapper::Config;
use CSG::Mapper::Exceptions;
use CSG::Storage::Slots;
use CSG::Types;

use Moose;
use Module::Load;

has '_conf' => (
  is      => 'ro',
  isa     => 'CSG::Mapper::Config',
  lazy    => 1,
  builder => '_build_conf',
);

has 'slot' => (
  is        => 'ro',
  isa       => 'Maybe[ValidSlotFSProject]',
  lazy      => 1,
  builder   => '_build_slot',
  predicate => 'has_slot',
  clearer   => 'clear_slot',
);

has 'cluster' => (
  is       => 'ro',
  isa      => 'ValidCluster',
  required => 1,
);

has 'build' => (
  is       => 'ro',
  isa      => 'Int',
  required => 1,
);

has 'record' => (
  is       => 'ro',
  isa      => 'CSG::Mapper::DB::Schema::Result::Sample',
  required => 1,
);

has 'prefix' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_prefix',
);

has 'project' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_project',
);

has 'host' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_host',
);

has 'center' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_center',
);

has 'pi' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_pi',
);

has 'run_dir' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_run_dir',
);

has 'filename' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_filename',
);

has 'sample_id' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_sample_id',
);

has 'incoming_path' => (
  is        => 'ro',
  isa       => 'Str',
  lazy      => 1,
  builder   => '_build_incoming_path',
  predicate => 'has_incoming_path',
);

has 'result_path' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_result_path',
);

has 'build_str' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_build_str',
);

has 'log_dir' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_log_dir',
);

has 'cram' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_cram',
);

has 'crai' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_crai',
);

has 'state_dir' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_state_dir',
);

before [qw(incoming_path result_path)] => sub {
  shift->slot;
};

after 'slot' => sub {
  my ($self) = @_;

  unless (defined $self->{slot}) {
    $self->clear_slot;
  }
};

sub _build_conf {
  return CSG::Mapper::Config->new(project => 'topmed');
}

sub _build_slot {
  my ($self) = @_;

  my $project = ucfirst(lc($self->project));
  my $class   = "CSG::Storage::SlotFS::$project";
  my $slot    = undef;

  try {
    load $class;

    # XXX - this is not right
    $slot = $class->find(
      name   => $self->sample_id,
      prefix => $self->prefix,
    );
  }
  catch {
    CSG::Mapper::Exceptions::Sample::SlotFailed->throw(error => $_);
  };

  return $slot;
}

sub _build_prefix {
  my ($self) = @_;
  return $self->_conf->get($self->cluster, 'prefix');
}

sub _build_project {
  return shift->record->project->name;
}

sub _build_host {
  return shift->record->host->name;
}

sub _build_center {
  return shift->record->center->name;
}

sub _build_pi {
  return shift->record->pi->name;
}

sub _build_run_dir {
  return shift->record->run_dir;
}

sub _build_filename {
  return shift->record->filename;
}

sub _build_sample_id {
  return shift->record->sample_id;
}

sub _build_incoming_path {
  my ($self) = @_;

  my ($bam, $cram) = (undef, undef);

  if ($self->has_slot) {
    $bam  = File::Spec->join($self->slot->incoming_path, $self->filename);
    $cram = File::Spec->join($self->slot->incoming_path, $self->sample_id . '.src.cram');
  } else {
    my $incoming_dir = $self->_conf->get($self->project, 'incoming_dir');
    my $backup_dir   = $self->_conf->get($self->project, 'backup_dir');
    my $base_dir = File::Spec->join($self->prefix, $self->host);

    # Original BAM path:
    # /<prefix>/<host>/<project_incoming_dir>/<center>/<run_dir>/<filename>
    $bam = File::Spec->join($base_dir, $incoming_dir, $self->center, $self->run_dir, $self->filename);

    # Backed up/squeezed path:
    # /<prefix>/<host>/<project_backup_dir>/<project_incoming_dir>/<center>/<run_dir>/<sample_id>.src.cram
    $cram = File::Spec->join($base_dir, $backup_dir, $incoming_dir, $self->center, $self->run_dir, $self->sample_id . '.src.cram');
  }

  return $bam  if -e $bam;
  return $cram if -e $cram;

  CSG::Mapper::Exceptions::Sample::NotFound->throw(bam_path => $bam, cram_path => $cram);
}

sub _build_result_path {
  my ($self) = @_;

  return File::Spec->join($self->slot->results_path, $self->build_str) if $self->has_slot;

  my $results_dir = $self->_conf->get($self->project, 'results_dir');
  return File::Spec->join($self->prefix, $self->host, $results_dir, $self->center, $self->pi, $self->build_str, $self->sample_id);
}

sub _build_cram {
  my ($self) = @_;
  return File::Spec->join($self->result_path, 'bams', $self->sample_id . '.recal.cram');
}

sub _build_crai {
  return shift->cram . '.crai';
}

sub _build_build_str {
  return 'hg' . shift->build;
}

sub _build_log_dir {
  my ($self) = @_;

  return $self->slot->log_path if $self->has_slot;

  my $log_dir = $self->_conf->get($self->project, 'log_dir');
  my $workdir = $self->_conf->get($self->project, 'workdir');

  return File::Spec->join($self->prefix, $workdir, $log_dir, $self->sample_id);
}

sub _build_state_dir {
  my ($self) = @_;

  return $self->slot->run_path if $self->has_slot;

  my $workdir = $self->_conf->get($self->project, 'workdir');
  my $run_dir = $self->_conf->get($self->project, 'run_dir');

  return File::Spec->join($self->prefix, $workdir, $run_dir, $self->sample_id);
}

sub is_complete {
  my ($self) = @_;

  return unless -e $self->cram;
  return if -z $self->cram;
  return unless -e $self->crai;

  my $cram_stat = File::Stat->new($self->cram);
  my $crai_stat = File::Stat->new($self->crai);

  return unless $crai_stat->mtime > $cram_stat->mtime;

  return 1;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
