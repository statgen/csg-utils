use utf8;
package CSG::Mapper::DB::Schema::Result::Result;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Mapper::DB::Schema::Result::Result

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

=head1 TABLE: C<results>

=cut

__PACKAGE__->table("results");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 sample_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 build

  data_type: 'varchar'
  default_value: 38
  is_nullable: 0
  size: 45

=head2 exported_at

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
  "sample_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "build",
  { data_type => "varchar", default_value => 38, is_nullable => 0, size => 45 },
  "exported_at",
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

=head1 UNIQUE CONSTRAINTS

=head2 C<index3>

=over 4

=item * L</sample_id>

=item * L</build>

=back

=cut

__PACKAGE__->add_unique_constraint("index3", ["sample_id", "build"]);

=head1 RELATIONS

=head2 jobs

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::Job>

=cut

__PACKAGE__->has_many(
  "jobs",
  "CSG::Mapper::DB::Schema::Result::Job",
  { "foreign.result_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 results_states_steps

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::ResultsStatesStep>

=cut

__PACKAGE__->has_many(
  "results_states_steps",
  "CSG::Mapper::DB::Schema::Result::ResultsStatesStep",
  { "foreign.result_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 sample

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Sample>

=cut

__PACKAGE__->belongs_to(
  "sample",
  "CSG::Mapper::DB::Schema::Result::Sample",
  { id => "sample_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2016-09-21 08:07:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:phDSY0AeCvoFSVP2FxNyjA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
#
sub status_line {
  my ($self) = @_;

  return sprintf
    q{ID: %-8s center: %-10s study: %-10s PI: %-15s Status: %-10s },
    $self->sample->sample_id,
    $self->sample->center->name,
    $self->sample->study->name,
    $self->sample->pi->name,
    $self->state->name;
}

sub cancel {
  my ($self) = @_;
  my $state = $self->result_source->schema->resultset('State')->find({name => 'cancelled'});
  $self->update({ state_id => $state->id });
  return;
}

sub current_state {
  my ($self) = @_;

  my $inside_rs = $self->results_states_steps->search();
  return $self->results_states_steps->find(
    {
      id => {'=' => $inside_rs->get_column('id')->max()},
    }
  )->state;
}

1;
