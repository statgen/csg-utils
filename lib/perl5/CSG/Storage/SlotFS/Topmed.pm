package CSG::Storage::SlotFS::Topmed;

use autodie qw(:all);
use Moose;

use Modern::Perl;
use Readonly;
use File::Path qw(make_path);

our $VERSION = '0.1';

Readonly::Scalar my $PROJECT => 'topmed';
Readonly::Array my @PATHS    => (qw(incoming results logs run info));

has 'size' => (
  is      => 'rw',
  isa     => 'Int',
  default => sub {
    return 0;
  },
);

has 'project' => (
  is      => 'ro',
  isa     => 'Str',
  default => sub {
    return $PROJECT;
  },
);

for my $path (@PATHS) {
  has "${path}_path" => (
    is      => 'ro',
    isa     => 'Str',
    lazy    => 1,
    builder => "_build_${path}_path",
  );

  eval "sub _build_${path}_path { return File::Spec->join(shift->path, $path); }";
}

sub initialize {
  my ($self) = @_;

  my @skel_dirs = map {File::Spec->join($self->path, $_)} @PATHS;

  make_path(@skel_dirs, {error => \my $err});

  if (@{$err}) {
    my $errstr;
    for (@{$err}) {
      my ($key, $value) = %{$_};
      $errstr = $value if $key eq '';
    }

    CSG::Storage::Slots::Exceptions::Sample::FailedSkeletonDirectory->throw(error => $errstr);
  }

  return;
}

sub to_string {
  return shift->path;
}

with qw(
  CSG::Storage::SlotFS::Roles::Sample
  CSG::Storage::SlotFS::Roles::Slot
  );

no Moose;
__PACKAGE__->meta->make_immutable;

1;

__END__


=head1 NAME

CSG::Storage::SlotFS::Topmed

=head1 VERSION

This documentation refers to CSG::Storage::SlotFS::Topmed version 0.1

=head1 SYNOPSIS

    use CSG::Storage::SlotFS::Topmed;

    my $topmed = CSG::Storage::SlotFS::Topmed->new(
      filename => 'foo.bam',
      name     => 'NWD123456',
    );

    $topmed->initialize(); # builds the topmed sample directory structure

    say $topmed->incoming_path; # path to the incoming directory for the sample

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
