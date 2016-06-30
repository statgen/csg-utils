use utf8;
package CSG::Storage::Slots::DB::Schema::Result::Slot;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Storage::Slots::DB::Schema::Result::Slot

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

=head1 TABLE: C<slots>

=cut

__PACKAGE__->table("slots");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 pool_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 sha1

  data_type: 'char'
  is_nullable: 1
  size: 40

=head2 size

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
  "pool_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "sha1",
  { data_type => "char", is_nullable => 1, size => 40 },
  "size",
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

=head1 UNIQUE CONSTRAINTS

=head2 C<index3>

=over 4

=item * L</name>

=item * L</pool_id>

=back

=cut

__PACKAGE__->add_unique_constraint("index3", ["name", "pool_id"]);

=head1 RELATIONS

=head2 pool

Type: belongs_to

Related object: L<CSG::Storage::Slots::DB::Schema::Result::Pool>

=cut

__PACKAGE__->belongs_to(
  "pool",
  "CSG::Storage::Slots::DB::Schema::Result::Pool",
  { id => "pool_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07043 @ 2016-03-23 08:59:19
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:6e1M73ACXFF6ISK5gc28QQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
use Number::Bytes::Human qw(format_bytes);

sub to_string {
  my ($self) = @_;
  return sprintf 'Id: %-5d Name: %-10s Project: %-10s Pool: %-10s Size: %s',
    $self->id,
    $self->name,
    $self->pool->project->name,
    $self->pool->name,
    format_bytes($self->size);
}

sub to_hashref {
  my ($self) = @_;
  return {
    id      => $self->id,
    name    => $self->name,
    project => $self->pool->project->name,
    pool    => $self->pool->name,
    size    => format_bytes($self->size),
  };
}

1;
