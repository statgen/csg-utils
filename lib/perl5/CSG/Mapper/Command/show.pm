package CSG::Mapper::Command::show;

use CSG::Mapper -command;
use CSG::Base qw(formats);
use CSG::Constants;
use CSG::Mapper::DB;

sub opt_spec {
  return (
    ['info',      'display basic job info'],
    ['meta-id=i', 'job meta id'],
    ['format=s',  'output format (yaml|txt)'],
    ['state|s=s', 'display results for a given state [submitted|requested|failed|cancelled|completed]'],
    ['build|b=s', 'reference build to show results for'],
  );
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  if ($opts->{info}) {
    unless ($opts->{meta_id}) {
      $self->usage_error('meta-id is required');
    }

    if ($opts->{format}) {
      unless ($opts->{format} =~ /yaml|txt/) {
        $self->usage_error('invalid output format');
      }
    }
  }

  if ($opts->{state}) {
    my $schema = CSG::Mapper::DB->new();

    unless ($schema->resultset('State')->find({name => $opts->{state}})) {
      $self->usage_error('invalid state');
    }

    unless ($opts->{build}) {
      $self->usage_error('build is required');
    }

    unless ($opts->{build} =~ /37|38/) {
      $self->usage_error('invalid build');
    }
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $schema = CSG::Mapper::DB->new();

  if ($opts->{info}) {
    my $format = $opts->{format} // 'yaml';
    my $meta   = $schema->resultset('Job')->find($opts->{meta_id});

    return $self->_info($meta, $format);
  }

  if ($opts->{state}) {
    my $results = $schema->resultset('Result')->search(
      {
        'me.build'   => $opts->{build},
        'state.name' => $opts->{state},
      },
      {
        join => 'state',
      }
    );

    for my $result ($results->all()) {
      say $result->status_line();
    }
  }
}

sub _info {
  my ($self, $meta, $format) = @_;

  my $info = {
    sample => {
      id        => $meta->result->sample->id,
      sample_id => $meta->result->sample->sample_id,
      center    => $meta->result->sample->center->name,
      study     => $meta->result->sample->study->name,
      pi        => $meta->result->sample->pi->name,
      host      => $meta->result->sample->host->name,
      filename  => $meta->result->sample->filename,
      run_dir   => $meta->result->sample->run_dir,
      state     => $meta->result->state->name,
      build     => $meta->result->build,
      fullpath  => $meta->result->sample->fullpath,
    },
    job => {
      id        => $meta->id,
      job_id    => $meta->job_id,
      cluster   => $meta->cluster,
      procs     => $meta->procs,
      memory    => $meta->memory,
      walltime  => $meta->walltime,
      node      => $meta->node,
      delay     => $meta->delay,
      submitted => ($meta->submitted_at) ? $meta->submitted_at->ymd . $SPACE . $meta->submitted_at->hms : $EMPTY,
      created   => $meta->created_at->ymd . $SPACE . $meta->created_at->hms,
    }
  };

  if ($format eq 'txt') {
    print Dumper $info;
  } else {
    print Dump($info);
  }

  return;
}

1;

__END__

=head1

CSG::Mapper::Command::show - show remapping jobs
