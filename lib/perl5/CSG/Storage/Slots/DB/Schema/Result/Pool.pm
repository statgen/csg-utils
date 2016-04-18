use utf8;
package CSG::Storage::Slots::DB::Schema::Result::Pool;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Storage::Slots::DB::Schema::Result::Pool

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

=head1 TABLE: C<pools>

=cut

__PACKAGE__->table("pools");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 project_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 type_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 hostname

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 size_used

  data_type: 'bigint'
  default_value: 0
  is_nullable: 0

=head2 size_total

  data_type: 'bigint'
  default_value: 0
  is_nullable: 0

=head2 path

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 active

  data_type: 'tinyint'
  default_value: 1
  is_nullable: 0

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
  "project_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "type_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "hostname",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "size_used",
  { data_type => "bigint", default_value => 0, is_nullable => 0 },
  "size_total",
  { data_type => "bigint", default_value => 0, is_nullable => 0 },
  "path",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "active",
  { data_type => "tinyint", default_value => 1, is_nullable => 0 },
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

=head2 C<index6>

=over 4

=item * L</name>

=item * L</hostname>

=back

=cut

__PACKAGE__->add_unique_constraint("index6", ["name", "hostname"]);

=head1 RELATIONS

=head2 project

Type: belongs_to

Related object: L<CSG::Storage::Slots::DB::Schema::Result::Project>

=cut

__PACKAGE__->belongs_to(
  "project",
  "CSG::Storage::Slots::DB::Schema::Result::Project",
  { id => "project_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 slots

Type: has_many

Related object: L<CSG::Storage::Slots::DB::Schema::Result::Slot>

=cut

__PACKAGE__->has_many(
  "slots",
  "CSG::Storage::Slots::DB::Schema::Result::Slot",
  { "foreign.pool_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 type

Type: belongs_to

Related object: L<CSG::Storage::Slots::DB::Schema::Result::Type>

=cut

__PACKAGE__->belongs_to(
  "type",
  "CSG::Storage::Slots::DB::Schema::Result::Type",
  { id => "type_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07043 @ 2016-03-31 14:00:57
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:NK7fAsuj7rlpX1cSCLQxtw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
sub max_available {
  my ($self) = @_;
  return $self->size_total - ($self->size_total * .10);
}

sub allocated {
  my ($self) = @_;
  my $allocated = $self->slots->search({},{+select => {sum => 'size', -as => 'allocated'}})->single;
  return $allocated->get_column('allocated') // 0;
}

sub is_full {
  my ($self) = @_;
  return $self->size_used >= $self->size_total;
}

sub over_threshold {
  my ($self, $size) = @_;
  return ($self->allocated + $size) > $self->max_available;
}

sub is_available {
  my ($self, $size) = @_;
  return unless $self->active;
  return if $self->is_full;
  return if $self->over_threshold($size);
  return 1;
}

1;
