package CSG::Mapper::Exceptions;

use Exception::Class (
  __PACKAGE__ . '::Job::BatchFileNotFound' => {
    description => 'Unable to locate batch file',
  },
  __PACKAGE__ . '::Job::BatchFileNotReadable' => {
    description => 'Unable to read batch file',
  },
  __PACKAGE__ . '::Job::SubmissionFailure' => {
    description => 'Failed to submit job',
  },
  __PACKAGE__ . '::Job::ProcessOutput' => {
    description => 'Failed to parse the job submission output',
    fields      => [qw(output)],
  },
  __PACKAGE__ . '::Job::NoJobToCancel' => {
    description => 'Can not cancel a job without a job id',
  },
  __PACKAGE__ . '::Job::CancellationFailure' => {
    description => 'Failed to cancel job',
  },
  __PACKAGE__ . '::Sample::NotFound' => {
    description => 'Sample not found on disk',
    fields      => [qw(bam_path cram_path)],
  },
  __PACKAGE__ . '::Sample::SlotFailed' => {
    description => 'Failed to load the slot for this sample',
  },
  __PACAKGE__ . '::Sample::FastqMismatch' => {
    description => 'Sample ID from list file does not match sample_id',
  },
  __PACKAGE__ . '::Sample::FastqNotFound' => {
    description => 'Fastq from list file does not exist on disk',
  },
  __PACKAGE__ . '::Sample::Fastq::MissingRG' => {
    description => 'Fastq list line read group field does not begin with @RG',
  },
  __PACKAGE__ . '::Sample::Fastq::MissingHeader' => {
    description => 'Fastq list line is missing a header field',
    fields      => [qw(header)],
  }

);

1;
