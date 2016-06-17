package Test::CSG::Storage::Slots;

use base qw(Test::CSG::Storage::Base);

use CSG::Base qw(test file);
use CSG::Storage::Slots;

use Digest::SHA qw(sha1_hex);
use Number::Bytes::Human qw(parse_bytes);

sub class {
  return 'CSG::Storage::Slots';
}

sub slot_prefix {
  return '/tmp/slots';
}

sub startup : Test(startup) {
  my ($self) = @_;
  make_path($self->slot_prefix) unless -e $self->slot_prefix;
  diag('Using prefix ' . $self->prefix);
}

sub setup : Test(setup) {
  my ($self) = @_;

  $ENV{SLOTS_DB} //= 'slots_test';

  my $config = CSG::Storage::Config->new();
  my $schema = CSG::Storage::Slots::DB->new();

  diag('Deploying schema to ' . $config->db);
  $schema->deploy({add_drop_table => 1});

  my $pools = YAML::LoadFile(File::Spec->join($self->fixture_path, 'pools.yml'));

  for my $pool (@{$pools}) {
    diag("Creating type: $pool->{type}");
    my $type = $schema->resultset('Type')->find_or_create({name => $pool->{type}});

    diag("Creating project: $pool->{project}");
    my $project = $schema->resultset('Project')->find_or_create({name => $pool->{project}});

    diag("Creating pool: $pool->{name}");
    my $pool = $schema->resultset('Pool')->find_or_create(
      {
        name       => $pool->{name},
        hostname   => $pool->{hostname},
        size_used  => parse_bytes($pool->{size_used}),
        size_total => parse_bytes($pool->{size_total}),
        path       => $pool->{path},
        type_id    => $type->id,
        project_id => $project->id,
      }
    );

    for my $slot (@{$pool->{slots}}) {
      diag("Creating slot: $slot->{name}");
      $pool->add_to_slots(
        {
          name => $slot->{name},
          size => $slot->{size},
        }
      );
    }
  }
}


sub test_new : Test(6) {
  my ($self) = @_;

  throws_ok {$self->class->new(prefix => $self->slot_prefix)} 'Moose::Exception::AttributeIsRequired', 'missing all params for new()';
  throws_ok {$self->class->new(prefix => $self->slot_prefix, name => 'foobar')} 'Moose::Exception::AttributeIsRequired', 'missing project and size for new()';
  throws_ok {$self->class->new(prefix => $self->slot_prefix, name => 'foobar', project => 'proj1')} 'Moose::Exception::AttributeIsRequired',
    'missing size for new()';

  lives_ok {$self->class->new(prefix => $self->slot_prefix, name => 'foobar', project => 'proj1', size => parse_bytes('100TB'))} 'all params given for new()';

  throws_ok {$self->class->new(prefix => $self->slot_prefix, name => 'foo', project => 'foo', size => parse_bytes('100TB'))}
  'CSG::Storage::Slots::Exceptions::Project::DoesNotExist', 'invalid project name';

  throws_ok {
    my $slot = $self->class->new(prefix => $self->slot_prefix, name => 'baz', project => 'proj1', size => parse_bytes('500TB'));
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
    prefix  => $self->slot_prefix,
  );

  like($slot->to_string, qr{^/tmp/slots/foo[1|2]\.localhost/working/slots[1|2]/8/8/4/3/foo}, 'to_string() path matches');
  like("$slot", qr{^/tmp/slots/foo[1|2]\.localhost/working/slots[1|2]/8/8/4/3/foo}, 'stringification path matches');
}

sub test_sha1 : Test(1) {
  my ($self) = @_;
  my $name   = 'foobarbaz';
  my $slot   = $self->class->new(prefix => $self->slot_prefix, name => $name, project => 'proj1', size => parse_bytes('200GB'));

  is($slot->sha1, sha1_hex($name), 'SHA1 should match');
}

sub test_parent : Test(4) {
  my ($self) = @_;
  lives_ok {$self->class->new(prefix => $self->slot_prefix, name => 'NA11931',      project => 'proj1', size => parse_bytes('100GB'))->path}
  'created parent slot';

  lives_ok {$self->class->new(prefix => $self->slot_prefix, name => 'NA11931-test', project => 'proj1', size => parse_bytes('100GB'))->path}
  'created child slot';

  throws_ok {$self->class->new(prefix => $self->slot_prefix, name => 'NA11931-foo',  project => 'proj1', size => parse_bytes('100GB'))->path}
  'CSG::Storage::Slots::Exceptions::Pools::NoPoolAvailable', 'not enough available space';

  throws_ok {$self->class->new(prefix => $self->slot_prefix, name => 'foobar-test', project => 'proj1', size => parse_bytes('100GB'))}
  'CSG::Storage::Slots::Exceptions::Slot::Parent::DoesNotExist', 'no parent slot found';
}

sub test_exists : Test(2) {
  my ($self) = @_;
  my $slot = $self->class->new(prefix => $self->slot_prefix, name => 'NA11931', project => 'proj1', size => parse_bytes('100GB'));
  $slot->path;
  ok($self->class->exists(name => 'NA11931', project => 'proj1'), 'Slot exists');
  ok(!$self->class->exists(name => 'jabberwocky', project => 'proj1'), 'Slot does not exist');
}

1;
