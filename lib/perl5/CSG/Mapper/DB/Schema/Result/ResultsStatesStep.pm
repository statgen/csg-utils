use utf8;
package CSG::Mapper::DB::Schema::Result::ResultsStatesStep;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Mapper::DB::Schema::Result::ResultsStatesStep

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

=head1 TABLE: C<results_states_steps>

=cut

__PACKAGE__->table("results_states_steps");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 result_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 state_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 step_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 job_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 created_at

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
  "state_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "step_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "job_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "created_at",
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

=head2 job

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Job>

=cut

__PACKAGE__->belongs_to(
  "job",
  "CSG::Mapper::DB::Schema::Result::Job",
  { id => "job_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
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

=head2 state

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::State>

=cut

__PACKAGE__->belongs_to(
  "state",
  "CSG::Mapper::DB::Schema::Result::State",
  { id => "state_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 step

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Step>

=cut

__PACKAGE__->belongs_to(
  "step",
  "CSG::Mapper::DB::Schema::Result::Step",
  { id => "step_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2016-10-10 14:35:38
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Q+J029beovb/s+GOqgtO/A


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
