package CSG::Storage::Slots::DB::Schema::ResultSet::Slot;

use base qw(DBIx::Class::ResultSet);

sub find_slot {
  my ($self, $name, $project) = @_;

  return $self->search(
    {
      'me.name'      => $name,
      'project.name' => $project,
    }, {
      join => {pool => 'project'},
    }
  )->first;
}

1;
