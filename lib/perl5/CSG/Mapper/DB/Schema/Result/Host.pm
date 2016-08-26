use utf8;
package CSG::Mapper::DB::Schema::Result::Host;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

CSG::Mapper::DB::Schema::Result::Host

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

=head1 TABLE: C<hosts>

=cut

__PACKAGE__->table("hosts");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

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
  "name",
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

=head2 fastqs

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::Fastq>

=cut

__PACKAGE__->has_many(
  "fastqs",
  "CSG::Mapper::DB::Schema::Result::Fastq",
  { "foreign.host_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 samples

Type: has_many

Related object: L<CSG::Mapper::DB::Schema::Result::Sample>

=cut

__PACKAGE__->has_many(
  "samples",
  "CSG::Mapper::DB::Schema::Result::Sample",
  { "foreign.host_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07043 @ 2016-08-24 10:19:38
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:YJ0GosEw1pVTwY6nkhMuvg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
