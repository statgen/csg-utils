package CSG::Types;

use CSG::Base;
use CSG::Constants qw(:mapping);

use Moose::Util::TypeConstraints;

subtype 'ValidCluster',
  as 'Str',
  where { $_ =~ /$VALID_CLUSTER_REGEXPS/ },
  message { 'is not a valid cluster type' };

subtype 'FileOnDisk',
  as 'Str',
  where {-s $_},
  message {'invalid file'};

subtype 'ValidJobFactory',
  as 'Object',
  where {
    $_->isa('CSG::Mapper::Job::Factory::Implementation::csg')
      or
    $_->isa('CSG::Mapper::Job::Factory::Implementation::flux')
  },
  message { 'is not a valid cluster job factory implementation' };

subtype 'Directory',
  as 'Str',
  where {-d $_},
  message { 'directory does not exist' };

no Moose::Util::TypeConstraints;

1;
