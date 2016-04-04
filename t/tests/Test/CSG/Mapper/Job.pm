package Test::CSG::Mapper::Job;

use base qw(Test::Class);
use CSG::Base qw(test);

use CSG::Mapper::Job;

sub class {
  return q{CSG::Mapper::Job};
}

sub cluster {
  return $ENV{TEST_CLUSTER};
}

sub startup : Test(startup) {
  my ($test) = @_;

  $test->{fixtures}->{clusters} = {
    csg => {
      job_id  => 12520268,
      state   => 'completed',
      elapsed => DateTime::Duration->new(
        days    => 2,
        hours   => 22,
        minutes => 30,
        seconds => 45,
      ),
    },
    flux => {
      job_id  => 16700562,
      state   => 'not_running',
      elapsed => DateTime::Duration->new(
        minutes => 2,
        seconds => 19,
      ),
    }
  };
}

sub setup : Test(setup => 6) {
  my ($test) = @_;

  for my $key (keys %{$test->{fixtures}->{clusters}}) {
    my $cluster = $test->{fixtures}->{clusters}->{$key};
    can_ok($test->class, 'new');
    my $job = $test->class->new(cluster => $key, job_id => $cluster->{job_id});
    isa_ok($job, $test->class);
    can_ok($job, 'job_id');
    $test->{fixtures}->{jobs}->{$key} = $job;
  }
}

sub test_job_id : Test(1) {
  my ($test) = @_;
  return unless $test->cluster;

  my $fixture = $test->{fixtures}->{clusters}->{$test->cluster};
  my $job     = $test->{fixtures}->{jobs}->{$test->cluster};

  is($fixture->{job_id}, $job->job_id);
}

sub test_submit : Test(1) {
  my ($test) = @_;
  return unless $test->cluster;

  my $fixture = $test->{fixtures}->{$test->cluster};
  my $job     = $test->{fixtures}->{jobs}->{$test->cluster};

  can_ok($job, 'submit');
}

1;
