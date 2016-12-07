use utf8;
package CSG::Mapper::DB::Schema::Result::Step;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Mapper::DB::Schema::Result::Step

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

=head1 TABLE: C<steps>

=cut

__PACKAGE__->table("steps");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  default_value: 'all'
  is_nullable: 0
  size: 45

=head2 parent_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "name",
  {
    data_type => "varchar",
    default_value => "all",
    is_nullable => 0,
    size => 45,
  },
  "parent_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");

=head1 UNIQUE CONSTRAINTS

=head2 C<name_UNIQUE>

=over 4

=item * L</name>

=back

=cut

__PACKAGE__->add_unique_constraint("name_UNIQUE", ["name"]);

=head1 RELATIONS

=head2 parent

Type: belongs_to

Related object: L<CSG::Mapper::DB::Schema::Result::Step>

=cut

__PACKAGE__->belongs_to(
  "parent",
  "CSG::Mapper::DB::Schema::Result::Step",
  { id => "parent_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
);

=head2 results_states_steps

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::ResultsStatesStep>

=cut

__PACKAGE__->has_many(
  "results_states_steps",
  "CSG::Mapper::DB::Schema::Result::ResultsStatesStep",
  { "foreign.step_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 steps

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::Step>

=cut

__PACKAGE__->has_many(
  "steps",
  "CSG::Mapper::DB::Schema::Result::Step",
  { "foreign.parent_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2016-12-07 11:27:12
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:f/A2BMEjqLWTeWeQew/OFA


# You can replace this text with custom code or comments, and it will be preserved on regeneration

sub has_parent {
  return defined shift->parent;
}

1;
