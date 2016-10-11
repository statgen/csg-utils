use utf8;
package CSG::Mapper::DB::Schema::Result::Job;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Mapper::DB::Schema::Result::Job

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

=head1 TABLE: C<jobs>

=cut

__PACKAGE__->table("jobs");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 result_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 job_id

  data_type: 'integer'
  is_nullable: 0

=head2 cluster

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 procs

  data_type: 'integer'
  is_nullable: 0

=head2 memory

  data_type: 'integer'
  is_nullable: 0

=head2 walltime

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 exit_code

  data_type: 'integer'
  is_nullable: 1

=head2 elapsed

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 node

  data_type: 'varchar'
  is_nullable: 1
  size: 45

=head2 delay

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 submitted_at

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 started_at

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 ended_at

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

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
  "result_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "job_id",
  { data_type => "integer", is_nullable => 0 },
  "cluster",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "procs",
  { data_type => "integer", is_nullable => 0 },
  "memory",
  { data_type => "integer", is_nullable => 0 },
  "walltime",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "exit_code",
  { data_type => "integer", is_nullable => 1 },
  "elapsed",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "node",
  { data_type => "varchar", is_nullable => 1, size => 45 },
  "delay",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "submitted_at",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "started_at",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "ended_at",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
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

=head2 logs

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::Log>

=cut

__PACKAGE__->has_many(
  "logs",
  "CSG::Mapper::DB::Schema::Result::Log",
  { "foreign.job_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 result

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Result>

=cut

__PACKAGE__->belongs_to(
  "result",
  "CSG::Mapper::DB::Schema::Result::Result",
  { id => "result_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 results_states_steps

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::ResultsStatesStep>

=cut

__PACKAGE__->has_many(
  "results_states_steps",
  "CSG::Mapper::DB::Schema::Result::ResultsStatesStep",
  { "foreign.job_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2016-10-11 08:06:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Ai/OQPQGjJAfa7XvxP7G7w

sub cancel {
  my ($self) = @_;

  $self->update(
    {
      ended_at  => $self->result_source->schema->now(),
      exit_code => '-1',
    }
  );

  $self->result->cancel();
}

sub result {
  my ($self) = @_;
  return $self->result_source->schema->resultset('ResultsStatesStep')->find(
    {
      id => {'=' =>  $self->results_states_steps->search()->get_column('id')->max()},
    }
  )->result;
}

1;
