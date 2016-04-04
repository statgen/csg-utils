package CSG::Mapper::Job::Factory;

use MooseX::AbstractFactory;

implementation_does [
  qw(
    CSG::Mapper::Job::Factory::Implementation::Requires
    )
];

implementation_class_via sub {
  q{CSG::Mapper::Job::Factory::Implementation::} . shift;
};

1;
