use utf8;
package CSG::Mapper::DB::Schema::Result::Sample;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Mapper::DB::Schema::Result::Sample

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 COMPONENTS LOADED

=over 4

=item * L<DBIx::Class::InflateColumn::DateTime>

=item * L<DBIx::Class::CSG::CreatedAt>

=back

=cut

__PACKAGE__->load_components("InflateColumn::DateTime", "CSG::CreatedAt");

=head1 TABLE: C<samples>

=cut

__PACKAGE__->table("samples");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 sample_id

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 center_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 study_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 pi_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 host_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 project_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 filename

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 run_dir

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 fullpath

  data_type: 'text'
  is_nullable: 0

=head2 year

  data_type: 'integer'
  is_nullable: 1

=head2 flagstat

  data_type: 'bigint'
  is_nullable: 1

=head2 ref_build

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 created_at

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 0

=head2 modified_at

  data_type: 'timestamp'
  datetime_undef_if_invalid: 1
  default_value: current_timestamp
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "sample_id",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "center_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "study_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "pi_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "host_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "project_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "filename",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "run_dir",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "fullpath",
  { data_type => "text", is_nullable => 0 },
  "year",
  { data_type => "integer", is_nullable => 1 },
  "flagstat",
  { data_type => "bigint", is_nullable => 1 },
  "ref_build",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "created_at",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 0,
  },
  "modified_at",
  {
    data_type => "timestamp",
    datetime_undef_if_invalid => 1,
    default_value => \"current_timestamp",
    is_nullable => 0,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");

=head1 RELATIONS

=head2 center

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Center>

=cut

__PACKAGE__->belongs_to(
  "center",
  "CSG::Mapper::DB::Schema::Result::Center",
  { id => "center_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 fastqs

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::Fastq>

=cut

__PACKAGE__->has_many(
  "fastqs",
  "CSG::Mapper::DB::Schema::Result::Fastq",
  { "foreign.sample_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 host

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Host>

=cut

__PACKAGE__->belongs_to(
  "host",
  "CSG::Mapper::DB::Schema::Result::Host",
  { id => "host_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 pi

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Pi>

=cut

__PACKAGE__->belongs_to(
  "pi",
  "CSG::Mapper::DB::Schema::Result::Pi",
  { id => "pi_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 project

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Project>

=cut

__PACKAGE__->belongs_to(
  "project",
  "CSG::Mapper::DB::Schema::Result::Project",
  { id => "project_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 results

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::Result>

=cut

__PACKAGE__->has_many(
  "results",
  "CSG::Mapper::DB::Schema::Result::Result",
  { "foreign.sample_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 study

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Study>

=cut

__PACKAGE__->belongs_to(
  "study",
  "CSG::Mapper::DB::Schema::Result::Study",
  { id => "study_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2016-11-07 08:16:17
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:cb3gwdeEmCyvoXPns/OoNg

use CSG::Constants;

sub has_fastqs {
  my ($self, $fastq) = @_;
  return $self->search_related('fastqs')->search({path => $fastq})->count;
}

sub has_fastqs_with_unpaired_reads {
  return shift->fastqs->search({path => {-not_like => '%interleaved%'}})->count;
}

sub has_unaligned_fastqs_in_read_group {
  my ($self, $rg) = @_;
  return $self->fastqs->search({read_group => $rg, aligned_at => undef})->count;
}

sub result_for_build {
  my ($self, $build) = @_;
  return $self->results->find({build => $build});
}

sub current_state {
  my ($self, $build) = @_;
}

sub is_requested {
  my ($self, $step, $build) = @_;
  my $result = $self->results->find({build => $build});

  return $FALSE unless $result;
  return $TRUE if $result->requested_step($step);
}

sub is_available {
  my ($self, $step, $build) = @_;

  my $result  = $self->results->find({build => $build});
  my $step_rs = $self->result_source->schema->resultset('Step')->find({name => $step});

  # XXX - has no result and step has no previous step
  return $TRUE if (not $result and not $step_rs->has_parent);

  # XXX - has no result and step has a previous step
  return $FALSE if (not $result and $step_rs->has_parent);

  # XXX - no results at all
  return $TRUE unless $result;

  # XXX - has to have completed the previous step
  return $FALSE unless $result->completed_previous_step($step);

  # XXX - has it completed the requested step already
  return $FALSE if $result->completed_step($step);

  # XXX - if the current state is requested then it needs processing
  return $TRUE if $result->requested_step($step);

  # XXX - states that mean manual intervention is required
  return $FALSE if $result->failed_step($step);
  return $FALSE if $result->cancelled_step($step);
  return $FALSE if $result->submitted_step($step);

  # XXX - passed everything else
  return $TRUE;
}

sub read_groups {
  return map {$_->read_group} shift->fastqs->search({}, {group_by => 'read_group'});
}

sub jobs_for_build {
  my ($self, $build) = @_;

  my $results = $self->result_for_build($build)->results_states_steps->search(
    {},
    {
      group_by => 'job_id',
      order_by => 'id desc'
    }
  );

  return $self->result_source->schema->resultset('Job')->search(
    {
      id => {'=' => [map {$_->job_id} $results->all]},
    },
    {
      order_by => 'id asc',
    }
  );
}

sub logs {
  my ($self, $build) = @_;

  my @logs = ();
  for my $job ($self->jobs_for_build($build)) {
    my $step   = $job->results_states_steps->first->step->name;
    my $job_id = $job->job_id;

    push @logs, map +{
      job_id    => $job_id,
      step      => $step,
      level     => uc($_->level),
      message   => $_->message,
      timestamp => $_->timestamp->strftime('%x %X'),
    }, $job->logs;
  }

  return wantarray ? @logs : \@logs;
}

1;
