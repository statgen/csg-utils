package CSG::Base;

use base qw(Import::Base);

our @IMPORT_MODULES = (
  qw(
    FindBin
    Modern::Perl
    Data::Dumper
    DateTime
    Readonly
    autodie
    Try::Tiny
    ),
  'English'         => [qw(-no_match_vars)],
  'List::MoreUtils' => [qw()],
  'Carp'            => [qw(croak)],
);

our %IMPORT_BUNDLES = (
  cmd => [
    'IPC::System::Simple' => [qw(run capture EXIT_ANY)],
    'System::Command',
  ],
  config => [
    qw(
      Config::Tiny
      )
  ],
  file => [
    qw(
      IO::All
      File::Spec
      File::Basename
      File::Stat
      Path::Class
      ),
    'File::Path'        => [qw(make_path remove_tree)],
    'File::Slurp::Tiny' => [qw(read_file read_lines)],
  ],
  formats => [
    qw(
      YAML
    )
  ],
  logging => [
    qw(
      Log::Dispatch
      Log::Dispatch::Screen
    )
  ],
  parsers => [
    qw(
      Class::CSV
    )
  ],
  templates => [
    qw(
      Template
    )
  ],
  test => [
    qw(
      Test::Most
      Test::More
      Test::Exception
      )
  ],
  www => [
    qw(
      URI
      URI::QueryParam
      JSON::MaybeXS
      Mojo::UserAgent
      )
  ]
);

1;
