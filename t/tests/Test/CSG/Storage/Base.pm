package Test::CSG::Storage::Base;

use base qw(Test::Class);
use Test::More;

use Modern::Perl;
use File::Path qw(remove_tree);
use File::Spec;
use File::Temp qw(tempdir);
use YAML qw(LoadFile);
use Number::Bytes::Human qw(parse_bytes);

use CSG::Storage::Config;
use CSG::Storage::Slots::DB;

my $PREFIX = tempdir();

sub fixture_path {
  return qq{$FindBin::Bin/../t/fixtures};
}

sub prefix {
  return $PREFIX;
}

sub _startup : Test(startup) {
  my ($self) = @_;

  my $config = CSG::Storage::Config->new();
  my $schema = CSG::Storage::Slots::DB->new();

  diag('Deploying schema to ' . $config->db);
  $schema->deploy({add_drop_table => 1});

  my $pools = LoadFile(File::Spec->join($self->fixture_path, 'pools.yml'));

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

sub _teardown : Test(teardown) {
  my ($self) = @_;
  diag('remove temporary directory ' . $self->prefix);
  remove_tree($self->prefix);
}

1;
