#!/bin/sh

export PERL_CARTON_PATH=[% settings.project_dir %]/local
export PERL5LIB=${PERL_CARTON_PATH}/lib/perl5:[% settings.project_dir %]/lib/perl5:${PERL5LIB}
export PATH=[% settings.project_dir %]/bin:${PERL_CARTON_PATH}/bin:${PATH}

META_ID=[% settings.meta_id %]
MAPPER_CMD=[% settings.mapper_cmd %]
MAPPER_LOG_CMD="$MAPPER_CMD log --meta-id $META_ID"
MAPPER_UPDATE_CMD="$MAPPER_CMD update --meta-id $META_ID"

$MAPPER_LOG_CMD --message 'sending fastq to the cloud'
$MAPPER_UPDATE_CMD --start --step cloud-align

# XXX - one idea, use gnu parallel to launch all fastqs at once
#
#       parallel -j <fastq_count> < fastq_align_commands.txt
#
# put this command for each fastq into a file:
#   cloud-align.sh [% fastq.read_group %] [% settings.out_dir %] [% fastq.path %]

cloud-align.sh [% fastq.read_group %] [% settings.out_dir %] [% fastq.path %]
if [ $? -ne 0 ]; then
  $MAPPER_LOG_CMD --message 'cloud alignment failed' --level critical
  $MAPPER_UPDATE_CMD --state failed
  exit 1
fi

CRAM="" # TODO - ?
if [ ! -e $CRAM ] ; then
  $MAPPER_LOG_CMD --message 'cloud alignment did not return a cram'
  $MAPPER_UPDATE_CMD --state failed
  exit 1
fi

# TODO - validate cram

$MAPPER_LOG_CMD --message 'removing fastq'
rm -v [% fastq.path %]

if [ $? -ne 0 ]; then
  $MAPPER_LOG_CMD --message 'failed to delete fastq'
  $MAPPER_UPDATE_CMD --state failed
  exit 1
fi

exit 0
