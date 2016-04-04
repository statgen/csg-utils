package CSG::Mapper::Command::freeze;

use CSG::Mapper -command;
use CSG::Base qw(parsers);
use CSG::Mapper::DB;
use CSG::Mapper::Sample;

sub opt_spec {
  return (['build|b=s', 'Reference build', {required => 1}]);
}

sub validate_args {
  my ($self, $opts, $args) = @_;

  unless ($self->app->global_options->{project}) {
    $self->usage_error('project name is required');
  }

  unless ($self->app->global_options->{cluster}) {
    $self->usage_error('cluster name is required');
  }
}

sub execute {
  my ($self, $opts, $args) = @_;

  my $fields  = [qw(sample_id center pi study cram)];
  my $csv     = Class::CSV->new(fields => $fields);
  my $schema  = CSG::Mapper::DB->new();
  my $project = $schema->resultset('Project')->find({name => $self->app->global_options->{project}});
  my @results = $project->completed_results($opts->{build})->all;

  $csv->add_line({map {$_ => $_} @{$fields}});


  for my $sample ($project->completed_results($opts->{build})->all) {
    my $sample_obj = CSG::Mapper::Sample->new(
      record  => $sample,
      cluster => $self->app->global_options->{cluster},
      build   => $opts->{build},
    );


    my $line_ref = {
      sample_id => $sample->sample_id,
      center    => $sample->center->name,
      pi        => $sample->pi->name,
      study     => $sample->study->name,
      cram      => $sample_obj->cram,
    };

    if ($opts->{build} eq '37') {
      # TODO - all paths for hg37 were built with previous toolsets and different paths
      $line_ref->{cram} = sprintf '/net/topmed/working/schelcj/results/%s/%s/%s/bams/%s.recal.cram',
        $sample->center->name,
        $sample->pi->name,
        $sample->sample_id,
        $sample->sample_id;

      next unless -e $line_ref->{cram};
      next unless -e $line_ref->{cram} . '.crai';

    } elsif ($opts->{build} eq '38') {
      next unless $sample_obj->is_complete;
    }

    $csv->add_line($line_ref);
  }

  $csv->print;
}

1;

__END__

=head1

CSG::Mapper::Command::freeze - cancel a remapping job
