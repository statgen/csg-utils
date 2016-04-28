requires 'Modern::Perl';
requires 'Readonly';
requires 'Try::Tiny';
requires 'IO::All';
requires 'Config::Tiny';
requires 'YAML';
requires 'Import::Base';
requires 'Template';
requires 'Class::CSV';
requires 'App::Cmd';
requires 'Exception::Class';
requires 'Module::Load';
requires 'Number::Bytes::Human';

requires 'Moose';
requires 'MooseX::AbstractFactory';

requires 'DateTime';
requires 'DateTime::Duration';
requires 'DateTime::Format::MySQL';

requires 'DBD::mysql';
requires 'DBD::SQLite';
requires 'DBIx::Class';
requires 'DBIx::Class::Schema::Loader';

requires 'Log::Dispatch';
requires 'Log::Dispatch::DBI';

requires 'System::Command';
requires 'IPC::System::Simple';

requires 'WWW::Mechanize';
requires 'Mojo::UserAgent';
requires 'JSON::MaybeXS';
requires 'IO::Socket::SSL';

requires 'File::Slurp::Tiny';
requires 'File::Stat';
requires 'File::Spec';
requires 'Filesys::DiskUsage';
requires 'Path::Class';

on 'test' => sub {
  requires 'SQL::Translator';
  requires 'Test::Class';
  requires 'Test::More';
  requires 'Test::Most';
  requires 'Test::Exception';
  requires 'Test::File';
};
