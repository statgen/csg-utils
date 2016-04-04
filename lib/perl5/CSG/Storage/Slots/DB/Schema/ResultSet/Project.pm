package CSG::Storage::Slots::DB::Schema::ResultSet::Project;

use base qw(DBIx::Class::ResultSet);

sub slots {
  my ($self, $name) = @_;

  return $self->search(
    {
      'me.name' => $name,
    },
    {
      join => {pools => 'slots'},
    }
  );
}

1;
