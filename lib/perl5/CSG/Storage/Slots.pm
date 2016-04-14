package CSG::Storage::Slots;

use Modern::Perl;
use Moose;
use Moose::Util::TypeConstraints;
use File::Spec;
use Digest::SHA qw(sha1_hex);
use overload '""' => sub {shift->to_string};

use CSG::Storage::Slots::DB;
use CSG::Storage::Slots::Exceptions;
use CSG::Types;

our $VERSION = "0.1";

subtype 'ValidProject',
  as 'Str',
  where {
    my $schema = CSG::Storage::Slots::DB->new();
    defined $schema->resultset('Project')->find({name => $_});
  },
  message {"Project name, $_, is not a valid project"};

has 'name'    => (is => 'ro', isa => 'Str',           required => 1);
has 'project' => (is => 'ro', isa => 'ValidProject',  required => 1);
has 'size'    => (is => 'rw', isa => 'ValidSlotSize', required => 1, trigger => \&_set_size);

has 'sha1' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_sha1');
has 'path' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_build_path');

has '_record' => (is => 'rw', isa => __PACKAGE__ . '::DB::Schema::Result::Slot', predicate => 'has_record');

around [qw(name project size path sha1)] => sub {
  my $orig = shift;
  my $self = shift;

  unless ($self->has_record) {
    my $schema  = CSG::Storage::Slots::DB->new();
    my $project = $schema->resultset('Project')->find({name => $self->{project}});
    my $pool    = $project->next_available_pool($self->{size});

    unless ($pool) {
      CSG::Storage::Slots::Exceptions::Pools::NoPoolAvailable->throw();
    }

    my $record  = $schema->resultset('Slot')->find_or_create(
      {
        name    => $self->{name},
        size    => $self->{size},
        sha1    => sha1_hex($self->{name}),
        pool_id => $pool->id,
      }
    );

    $self->_record($record);
  }

  return $self->$orig(@_);
};

sub _set_size {
  my ($self, $new, $old) = @_;

  if ($self->has_record) {
    $self->_record->update({size => $new});
  }
}

sub _build_sha1 {
  return shift->_record->sha1;
}

sub _build_path {
  my ($self) = @_;
  my $pool = $self->_record->pool;
  return File::Spec->join($pool->hostname, $pool->path, (split(//, $self->sha1))[0 .. 3], $self->name);
}

sub to_string {
  return shift->path;
}

sub find {
  my $class  = shift;
  my %params = @_;
  my $schema = CSG::Storage::Slots::DB->new();
  my $slot   = $schema->resultset('Slot')->find_slot($params{name}, $params{project});

  return unless $slot;

  return $class->new(
    name    => $slot->name,
    project => $slot->pool->project->name,
    size    => $slot->size,
    _record => $slot,
  );
}

sub find_or_create {
  my $class = shift;
  return $class->find(@_) // $class->new(@_);
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 NAME

 CSG::Storage::Slots

=head1 VERSION

This documentation refers to CSG::Storage::Slots version 0.1

=head1 SYNOPSIS

    use CSG::Storage::Slots;

    my $slot = CSG::Storage::Slots->new(name => 'foo', project => 'bar', size => '2000000000');

    say $slot->path;

=head1 DESCRIPTION



=head1 DIAGNOSTICS



=head1 CONFIGURATION AND ENVIRONMENT



=head1 DEPENDENCIES



=head1 INCOMPATIBILITIES



=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.
Please report problems to Chris Scheller <schelcj@umich.edu>
Patches are welcome.

=head1 AUTHOR

Chris Scheller <schelcj@umich.edu>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2016 Regents of the University of Michigan. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
