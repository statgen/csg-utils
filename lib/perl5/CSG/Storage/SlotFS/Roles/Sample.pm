package CSG::Storage::SlotFS::Roles::Sample;

use autodie qw(:all);
use Moose::Role;

use Modern::Perl;
use IPC::System::Simple qw(capture);

use CSG::Storage::Types;

has 'filename' => (
  is        => 'rw',
  isa       => 'Maybe[ValidFile]',
  predicate => 'has_filename'
);

has 'sample_id' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_sample_id'
);

sub _build_sample_id {
  my ($self) = @_;
  return unless $self->has_filename;
  return capture(sprintf q{samtools view -H %s|grep '^@RG'|grep -o 'SM:\S*'|sort -u|cut -d \: -f 2}, $self->filename);
}

1;
