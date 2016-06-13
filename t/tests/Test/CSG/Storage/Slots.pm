package Test::CSG::Storage::Slots;

use base qw(Test::CSG::Storage::Base);

use CSG::Base qw(test);
use CSG::Storage::Slots;

use Digest::SHA qw(sha1_hex);
use Number::Bytes::Human qw(parse_bytes);

sub class {
  return 'CSG::Storage::Slots';
}

sub test_new : Test(6) {
  my ($self) = @_;

  throws_ok {$self->class->new()} 'Moose::Exception::AttributeIsRequired', 'missing all params for new()';
  throws_ok {$self->class->new(name => 'foobar')} 'Moose::Exception::AttributeIsRequired', 'missing project and size for new()';
  throws_ok {$self->class->new(name => 'foobar', project => 'proj1')} 'Moose::Exception::AttributeIsRequired',
    'missing size for new()';

  lives_ok {$self->class->new(name => 'foobar', project => 'proj1', size => parse_bytes('100TB'))} 'all params given for new()';

  throws_ok {$self->class->new(name => 'foo', project => 'foo', size => parse_bytes('100TB'))}
  'CSG::Storage::Slots::Exceptions::Project::DoesNotExist', 'invalid project name';

  throws_ok {
    my $slot = $self->class->new(name => 'baz', project => 'proj1', size => parse_bytes('500TB'));
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
    prefix  => '/tmp',
  );

  like($slot->to_string, qr{^/tmp/foo[1|2]\.localhost/working/slots[1|2]/8/8/4/3/foo}, 'to_string() path matches');
  like("$slot", qr{^/tmp/foo[1|2]\.localhost/working/slots[1|2]/8/8/4/3/foo}, 'stringification path matches');
}

sub test_sha1 : Test(1) {
  my ($self) = @_;
  my $name   = 'foobarbaz';
  my $slot   = $self->class->new(name => $name, project => 'proj1', size => parse_bytes('200GB'));

  is($slot->sha1, sha1_hex($name), 'SHA1 should match');
}

1;
