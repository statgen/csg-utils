package CSG::Storage::Slots::Exceptions;

use Exception::Class (
  __PACKAGE__ . '::SlotExists' => {
    description => 'slot directory already exists',
  },
  __PACKAGE__ . '::Slot::DoesNotExist' => {
    description => 'unable to locate requested slot',
  },
  __PACKAGE__ . '::Slot::Parent::DoesNotExist' => {
    description => 'no parent slot found',
  },
  __PACKAGE__ . '::Sample::FailedCopy' => {
    description => 'failed to copy sample to incoming directory',
    fields      => [qw(error)],
  },
  __PACKAGE__ . '::Sample::FailedSkeletonDirectory' => {
    description => 'failed to create skeleton directory sturcture',
    fields      => [qw(error)],
  },
  __PACKAGE__ . '::Pools::NoPoolAvailable' => {
    description => 'unable to find pool with enough available space',
  },
  __PACKAGE__ . '::Project::DoesNotExist' => {
    description => 'project does not exist',
  },
);

1;
