package Test::CSG::Mapper::Test;

use base qw(Test::Class);
use CSG::Base qw(test);

use CSG::Mapper::Util qw(:all);

sub class {
  return 'CSG::Mapper::Util';
}

sub cluster {
  return $ENV{TEST_CLUSTER} // 'csg';
}

sub test_detect_cluster : Test(no_plan) {
  is(detect_cluster(), shift->cluster, 'found correct distro');
}

1;
