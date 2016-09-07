use utf8;
package CSG::Mapper::DB::Schema::Result::Fastq;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Mapper::DB::Schema::Result::Fastq

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

=head1 TABLE: C<fastqs>

=cut

__PACKAGE__->table("fastqs");

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

=head2 path

  data_type: 'text'
  is_nullable: 0

=head2 read_group

  data_type: 'text'
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
  "sample_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "build",
  { data_type => "varchar", default_value => 38, is_nullable => 0, size => 45 },
  "path",
  { data_type => "text", is_nullable => 0 },
  "read_group",
  { data_type => "text", is_nullable => 0 },
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


# Created by DBIx::Class::Schema::Loader v0.07043 @ 2016-09-07 08:42:41
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:vbRoLSsp1lQ8XlRbHhaFVw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
