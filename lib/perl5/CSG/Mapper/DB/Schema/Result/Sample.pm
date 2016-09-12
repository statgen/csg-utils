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


# Created by DBIx::Class::Schema::Loader v0.07043 @ 2016-08-24 10:19:38
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:DHrlmMLVhJxoRPq4xJdhow

sub has_fastqs {
  my ($self, $fastq) = @_;
  return $self->search_related('fastqs')->search({path => $fastq})->count;
}

1;
