#!/bin/bash
#
# Usage:
#   gce-align.sh [paired] [sample_id] [read_group_ident] [read_group] [output_dir] [input_files]
#
#   [paired]           - either 0 or 1 for unpaired or paired reads
#   [sample_id]        - sample identifier (e.g. nwd123456)
#   [read_group_ident] - identification string to describe the read group for vm creatation (e.g rg-0)
#   [read_group]       - read group string from the .list file during bam2fastq
#   [output_dir]       - directory, or google bucket, where cram files will be written
#   [input_files]      - full path to fastq file(s). can be absolute path or google bucket uri (e.g. gs://topmed-fastqs/)
#
PAIRED=$1
shift
SAMPLE_ID=$1
shift
RG_IDENT=$1
shift
RG_LINE=$1
shift
OUTPUT_DIR=$1
shift
INPUT_FILES=$@

set -o pipefail

case "$PAIRED" in
  0) BWA_PAIRED_OPT="" ;;
  1) BWA_PAIRED_OPT="-p" ;;
esac

MACHINE_NAME="align-${SAMPLE_ID}-${RG_IDENT}"
MACHINE_TYPE_OPTS="--custom-cpu 32 --custom-memory 64GiB"

gcloud compute instances create --scopes storage-full --image ubuntu-1404-alignment $MACHINE_TYPE_OPTS --boot-disk-size 300 --preemptible $MACHINE_NAME

EXIT_STATUS=$?
if [[ $EXIT_STATUS == 0 ]]
then
  # Give ssh daemon some time to get running.
  sleep 30s

  #COMMAND="set -o pipefail; bwa mem -t 32 -K 100000000 -Y -p -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samtools sort -@ 32 -m 2000M --reference /home/alignment/ref/hs38DH.fa -O cram -o /home/alignment/output.cram -T /home/alignment/sort.temp -"
  COMMAND="set -o pipefail; rm -f /home/alignment/output.cram /home/alignment/output.cram.ok; bwa mem -t 32 -K 100000000 -Y $BWA_PAIRED_OPT -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samblaster -a --addMateTags | samtools view -@ 32 -T /home/alignment/ref/hs38DH.fa -C -o /home/alignment/output.cram - && touch /home/alignment/output.cram.ok"

  CONTAINER_ID=$(gcloud compute ssh $MACHINE_NAME -- sudo docker create -v "/home/alignment:/home/alignment" statgen/alignment /bin/bash -c \""$COMMAND"\")

  EXIT_STATUS=$?
  if [[ $EXIT_STATUS == 0 ]]
  then
    for f in $INPUT_FILES
    do
      echo "[$(date)] Downloading $f"
      START_TIME=$(date +%s)
      gcloud compute ssh $MACHINE_NAME -- gsutil cp $f /home/alignment/input.fastq.gz
      EXIT_STATUS=$?
      echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"

      RETRY_COUNTER=0
      while [[ $EXIT_STATUS == 0 && $RETRY_COUNTER -lt 5 ]]
      do

        gcloud compute ssh $MACHINE_NAME -- sudo docker start $CONTAINER_ID

        CONTAINER_IS_RUNNING=1
        FAILED_CONTAINER_POLL_COUNT=0
        while [[ $CONTAINER_IS_RUNNING != 0 && $FAILED_CONTAINER_POLL_COUNT -lt 5 ]]
        do
          sleep 60s

          CONTAINER_STATUS=$(gcloud compute ssh $MACHINE_NAME -- sudo docker ps -al --format {{.Status}})
          if [[ $CONTAINER_STATUS =~ Exited\ \((.*)\) ]]
          then
            CONTAINER_IS_RUNNING=0
            EXIT_STATUS=${BASH_REMATCH[1]}

            echo "[$(date)] Fetching logs"
            gcloud compute ssh $MACHINE_NAME --command "sudo docker logs $CONTAINER_ID"
            if [[ $? != 0 ]]
            then
              echo "[$(date)] Fetching logs FAILED!"
            fi
          fi

          # if [[ $(gcloud compute instances list $MACHINE_NAME | grep RUNNING | wc -l) == 0 ]]
          # then
          #   MACHINE_IS_RUNNING=0
          # fi
          if [[ $CONTAINER_STATUS ]]
          then
            FAILED_CONTAINER_POLL_COUNT=0
          else
            let FAILED_CONTAINER_POLL_COUNT++
            echo "[$(date)] Failed Count: $FAILED_CONTAINER_POLL_COUNT"
          fi
        done

        if [[ $CONTAINER_IS_RUNNING == 0 && $EXIT_STATUS == 0 ]]
        then
          OUTPUT_FILE=$OUTPUT_DIR"/"$(basename $f .fastq.gz)".cram"
          echo "[$(date)] Uploading $OUTPUT_FILE"
          START_TIME=$(date +%s)
          gcloud compute ssh $MACHINE_NAME -- gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp /home/alignment/output.cram $OUTPUT_FILE && gcloud compute ssh $MACHINE_NAME -- gsutil cp /home/alignment/output.cram.ok $OUTPUT_FILE".ok"
          EXIT_STATUS=$?
          echo "[$(date)] Upload exit status: $EXIT_STATUS"
          echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"

          break
        elif [[ $FAILED_CONTAINER_POLL_COUNT == 5 && $RETRY_COUNTER -lt 4 ]]
        then
          #sleep 300s
          echo "[$(date)] Machine stopped. Restarting $MACHINE_NAME"
          gcloud compute instances start $MACHINE_NAME
        fi

        let RETRY_COUNTER++
      done
    done
  fi
fi

# Machine may exist even if "docker-machine create" fails.
echo "[$(date)] Cleaning up ..."
gcloud compute instances delete --quiet $MACHINE_NAME

exit $EXIT_STATUS
