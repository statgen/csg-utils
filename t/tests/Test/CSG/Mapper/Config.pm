package Test::CSG::Mapper::Test;

use base qw(Test::Class);
use CSG::Base qw(test);

use CSG::Mapper::Config;

sub class {
  return 'CSG::Mapper::Config';
}

sub startup : Test(startup => 2) {
  my ($test) = @_;

  my $fixture_path = qq{$FindBin::Bin/../t/fixtures/configs};

  $ENV{CSG_MAPPING_CONF} = qq{$FindBin::Bin/../config/mapper.ini};

  my $config = $test->class->new(
    project     => 'topmed',
    _config_dir => $fixture_path,
  );

  isa_ok($config, $test->class);
  is($config->_config_dir, $fixture_path, 'config directory path matches');
  $test->{config} = $config;
}

sub test_dsn : Test(1) {
  my ($test) = @_;
  my $config = $test->{config};

  is($config->dsn, 'dbi:mysql:database=baz;host=localhost;port=3306', 'dsn matches');
}

sub test_project : Test(2) {
  my ($test) = @_;
  dies_ok(sub {$test->class->new(project => 'foo')}, 'expected to die');
  lives_ok(sub {$test->class->new(project => 'topmed')}, 'valid project lives');
}

sub test_get : Test(8) {
  my ($test) = @_;
  my $config = $test->{config};

  is($config->get('topmed',    'procs'),         6,                        'procs matches');
  is($config->get('topmed',    'walltime'),      '672:00:00',              'walltime matches');
  is($config->get('topmed',    'build'),         38,                       'build matches');
  is($config->get('pipelines', 'uw'),            'cleanUpBam2fastq',       'pipeline matches');
  is($config->get('csg',       'gotcloud_conf'), 'config/gotcloud.conf.csg',  'gotcloud conf matches for cluster csg');
  is($config->get('flux',      'gotcloud_conf'), 'config/gotcloud.conf.flux', 'gotcloud conf matches for cluster flux');

  is($config->get('foo',    'bar'), undef, 'expected undef for non-existant category');
  is($config->get('topmed', 'bar'), undef, 'expected undef for non-existant value');
}

sub test_overrides : Test(no_plan) {
  my ($test) = @_;
  my $config = $test->{config};

  is($config->get('csg', 'align_procs'), 8, 'project specific cores for align step match');
}

1;
