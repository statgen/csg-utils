package Test::CSG::Storage::Slots;

use base qw(Test::CSG::Storage::Base);

use CSG::Base qw(test file);
use CSG::Storage::Slots;

use Digest::SHA qw(sha1_hex);
use Number::Bytes::Human qw(parse_bytes);

sub class {
  return 'CSG::Storage::Slots';
}

sub prefix {
  return '/tmp/slots';
}

sub setup_prefix : Test(setup) {
  my ($self) = @_;
  diag($self->prefix);
  make_path($self->prefix) unless -e $self->prefix;
}

sub test_new : Test(6) {
  my ($self) = @_;

  throws_ok {$self->class->new(prefix => $self->prefix)} 'Moose::Exception::AttributeIsRequired', 'missing all params for new()';
  throws_ok {$self->class->new(prefix => $self->prefix, name => 'foobar')} 'Moose::Exception::AttributeIsRequired', 'missing project and size for new()';
  throws_ok {$self->class->new(prefix => $self->prefix, name => 'foobar', project => 'proj1')} 'Moose::Exception::AttributeIsRequired',
    'missing size for new()';

  lives_ok {$self->class->new(prefix => $self->prefix, name => 'foobar', project => 'proj1', size => parse_bytes('100TB'))} 'all params given for new()';

  throws_ok {$self->class->new(prefix => $self->prefix, name => 'foo', project => 'foo', size => parse_bytes('100TB'))}
  'CSG::Storage::Slots::Exceptions::Project::DoesNotExist', 'invalid project name';

  throws_ok {
    my $slot = $self->class->new(prefix => $self->prefix, name => 'baz', project => 'proj1', size => parse_bytes('500TB'));
    $slot->size;
  }
  'CSG::Storage::Slots::Exceptions::Pools::NoPoolAvailable', 'not enough space in exising pools';
}

sub test_path : Test(2) {
  my ($self) = @_;

  my $slot = $self->class->new(
    name    => 'foobar',
    project => 'proj1',
    size    => parse_bytes('300GB'),
    prefix  => $self->prefix,
  );

  like($slot->to_string, qr{^/tmp/slots/foo[1|2]\.localhost/working/slots[1|2]/8/8/4/3/foo}, 'to_string() path matches');
  like("$slot", qr{^/tmp/slots/foo[1|2]\.localhost/working/slots[1|2]/8/8/4/3/foo}, 'stringification path matches');
}

sub test_sha1 : Test(1) {
  my ($self) = @_;
  my $name   = 'foobarbaz';
  my $slot   = $self->class->new(prefix => $self->prefix, name => $name, project => 'proj1', size => parse_bytes('200GB'));

  is($slot->sha1, sha1_hex($name), 'SHA1 should match');
}

sub test_parent : Test(4) {
  my ($self) = @_;
  lives_ok {$self->class->new(prefix => $self->prefix, name => 'NA11931',      project => 'proj1', size => parse_bytes('100GB'))->path}
  'created parent slot';

  lives_ok {$self->class->new(prefix => $self->prefix, name => 'NA11931-test', project => 'proj1', size => parse_bytes('100GB'))->path}
  'created child slot';

  throws_ok {$self->class->new(prefix => $self->prefix, name => 'NA11931-foo',  project => 'proj1', size => parse_bytes('100GB'))->path}
  'CSG::Storage::Slots::Exceptions::Pools::NoPoolAvailable', 'not enough available space';

  throws_ok {$self->class->new(prefix => $self->prefix, name => 'foobar-test', project => 'proj1', size => parse_bytes('100GB'))}
  'CSG::Storage::Slots::Exceptions::Slot::Parent::DoesNotExist', 'no parent slot found';
}

1;
