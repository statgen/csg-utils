#!/bin/bash
set -o pipefail

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
INPUT_FILES_COUNT=$#

case "$PAIRED" in
  0) BWA_PAIRED_OPT="" ;;
  1) BWA_PAIRED_OPT="-p" ;;
esac

ZONES[0]="us-central1-b"
ZONES[1]="us-central1-c"
ZONES[2]="us-central1-f"
MACHINE_ZONE=${ZONES[$[ $RANDOM % 3 ]]}

MACHINE_NAME=$(echo "align-${SAMPLE_ID}-${RG_IDENT}" | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)
MACHINE_TYPE_OPTS="--custom-cpu 32 --custom-memory 64GiB"

EXIT_STATUS=0
RETRY_COUNTER=0
while [[ $# -gt 0 && $EXIT_STATUS == 0 && $RETRY_COUNTER -lt 5 ]]
do

  gcloud compute instances create --zone $MACHINE_ZONE --scopes storage-full --image ubuntu-1404-alignment $MACHINE_TYPE_OPTS --boot-disk-size 300 --preemptible $MACHINE_NAME

  EXIT_STATUS=$?
  if [[ $EXIT_STATUS == 0 ]]
  then
    # Give ssh daemon some time to get running.
    sleep 30s

    #COMMAND="set -o pipefail; bwa mem -t 32 -K 100000000 -Y -p -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samtools sort -@ 32 -m 2000M --reference /home/alignment/ref/hs38DH.fa -O cram -o /home/alignment/output.cram -T /home/alignment/sort.temp -"
    COMMAND="set -o pipefail; rm -f /home/alignment/output.cram /home/alignment/output.cram.ok; bwa mem -t 32 -K 100000000 -Y $BWA_PAIRED_OPT -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samblaster -a --addMateTags | samtools view -@ 32 -T /home/alignment/ref/hs38DH.fa -C -o /home/alignment/output.cram - && touch /home/alignment/output.cram.ok"

    CONTAINER_ID=$(gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- sudo docker create -v "/home/alignment:/home/alignment" statgen/alignment /bin/bash -c \""$COMMAND"\")

    EXIT_STATUS=$?
    if [[ $EXIT_STATUS == 0 ]]
    then 
      while [[ $# -gt 0 ]]
      do
        F=$1
        echo "[$(date)] Downloading "$F" ..."
        START_TIME=$(date +%s)
        gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- gsutil cp $F /home/alignment/input.fastq.gz
        EXIT_STATUS=$?
        echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"
    
        gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- sudo docker start $CONTAINER_ID

        CONTAINER_IS_RUNNING=1
        FAILED_CONTAINER_POLL_COUNT=0
        while [[ $CONTAINER_IS_RUNNING != 0 && $FAILED_CONTAINER_POLL_COUNT -lt 5 ]]
        do
          sleep 60s

          CONTAINER_STATUS=$(gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- sudo docker ps -al --format {{.Status}})
          if [[ $CONTAINER_STATUS =~ Exited\ \((.*)\) ]]
          then
            CONTAINER_IS_RUNNING=0
            EXIT_STATUS=${BASH_REMATCH[1]}

            echo "[$(date)] Fetching logs ..."
            gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME --command "sudo docker logs $CONTAINER_ID"
            if [[ $? != 0 ]]
            then
              echo "[$(date)] Fetching logs FAILED!"
            fi
          fi

          if [[ $CONTAINER_STATUS ]]
          then
            FAILED_CONTAINER_POLL_COUNT=0
          else
            let FAILED_CONTAINER_POLL_COUNT++
            echo "[$(date)] Failed Count: "$FAILED_CONTAINER_POLL_COUNT
          fi
        done

        if [[ $CONTAINER_IS_RUNNING == 0 && $EXIT_STATUS == 0 ]]
        then
          OUTPUT_FILE=$OUTPUT_DIR"/"$(basename $F .fastq.gz)".cram"
          echo "[$(date)] Uploading "$OUTPUT_FILE" ..."
          START_TIME=$(date +%s)
          gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp /home/alignment/output.cram $OUTPUT_FILE && gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- gsutil cp /home/alignment/output.cram.ok $OUTPUT_FILE".ok"
          EXIT_STATUS=$?
          echo "[$(date)] Upload exit status: "$EXIT_STATUS
          echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"
          
          if [[ $RETRY_COUNTER -gt 0 ]]
          then
            let RETRY_COUNTER--
          fi

          shift #pops file off of input list
        elif [[ $FAILED_CONTAINER_POLL_COUNT == 5 ]]
        then
          echo "[$(date)] Machine stopped: "$MACHINE_NAME
          break
        fi
      done
    fi
  fi

  # Machine may exist even if "docker-machine create" fails.
  echo "[$(date)] Cleaning up ..."
  gcloud compute instances delete --zone $MACHINE_ZONE --quiet $MACHINE_NAME
  
  if [[ $# -gt 0 ]]
  then
    sleep $(( $RETRY_COUNTER * 600 ))"s"
  fi

  let RETRY_COUNTER++
done

if [[ $# -gt 0 ]]
then
  echo "[$(date)] Giving up at "$1
fi

echo "Done."

exit $EXIT_STATUS