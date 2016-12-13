#!/bin/bash
set -o pipefail

ZONES=("us-central1-a" "us-central1-b" "us-central1-c" "us-central1-f")
MACHINE_ZONE=${ZONES[$[ $RANDOM % ${#ZONES[@]} ]]}

MACHINE_NAME=$(echo "pre-align-"$(basename $1 | cut -f 1 -d '.') | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)
#MACHINE_TYPE_OPTS="--machine-type n1-standard-2"
MACHINE_TYPE_OPTS="--custom-cpu 1 --custom-memory 6656MiB"

gcloud compute instances create --zone $MACHINE_ZONE --scopes storage-full --image ubuntu-1404-alignment $MACHINE_TYPE_OPTS --boot-disk-size 300 $MACHINE_NAME

EXIT_STATUS=$?
if [[ $EXIT_STATUS == 0 ]]
then
  # Give ssh daemon some time to get running.
  sleep 30s  

  echo "[$(date)] Downloading input crams"
  START_TIME=$(date +%s)
  INPUT_FILE="/home/alignment/"$(basename $1)
  OUT_BASE="/home/alignment/"$(basename $1 | cut -f 1 -d '.')
 
  gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- gsutil cp $1 $INPUT_FILE
  EXIT_STATUS=$?
  echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"

  if [[ $EXIT_STATUS == 0 ]]
  then
    REMOTE_COMMAND="set -euo pipefail &&
      export REF_PATH=/home/alignment/ref/md5/%2s/%2s/%s &&
      samtools view -uh -F 0x900 $INPUT_FILE \
      | bam-ext-mem-sort-manager squeeze --in -.ubam --keepDups --rmTags AS:i,BD:Z,BI:Z,XS:i,MC:Z,MD:Z,NM:i,MQ:i --out -.ubam \
      | samtools sort -l 1 -@ 1 -m 4000M -n -T /home/alignment/sort_tmp - \
      | samtools fixmate - - \
      | bam-ext-mem-sort-manager bam2fastq --in -.bam --outBase $OUT_BASE --maxRecordLimitPerFq 20000000 --sortByReadNameOnTheFly --readname --gzip
      fastq_reads=\\\$((\\\$(zcat /home/alignment/*.fastq.gz | wc -l) / 4))
      samtools flagstat $INPUT_FILE > /home/alignment/cram_flagstat.txt
      cram_reads=\\\$(grep 'paired in sequencing' /home/alignment/cram_flagstat.txt | awk {'print \\\$1'})
      echo 'In read count: '\\\$cram_reads
      echo 'Out read count: '\\\$fastq_reads
      if [[ \\\$fastq_reads != \\\$cram_reads ]]; then echo FAILED; exit -1; fi"

    echo "[$(date)] Creating container"
    CONTAINER_ID=$(gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- sudo docker create -v "/home/alignment:/home/alignment" statgen/alignment /bin/bash -c \""$REMOTE_COMMAND"\")

    echo "[$(date)] Starting container"
    gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- sudo docker start $CONTAINER_ID

    CONTAINER_IS_RUNNING=1
    FAILED_CONTAINER_POLL_COUNT=0
    while [[ $CONTAINER_IS_RUNNING != 0 && $FAILED_CONTAINER_POLL_COUNT -lt 5 ]]
    do
      sleep 180s

      CONTAINER_STATUS=$(gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- sudo docker ps -al --format {{.Status}})
      if [[ $CONTAINER_STATUS =~ Exited\ \((.*)\) ]]
      then
        CONTAINER_IS_RUNNING=0
        EXIT_STATUS=${BASH_REMATCH[1]}

        echo "[$(date)] Fetching logs"
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
        echo "[$(date)] Failed Count: $FAILED_CONTAINER_POLL_COUNT"
      fi
    done

    if [[ $CONTAINER_IS_RUNNING == 0 && $EXIT_STATUS == 0 ]]
    then
      OUTPUT_DIR=$2"/"
      echo "[$(date)] Uploading $OUTPUT_FILE"
      START_TIME=$(date +%s)
      gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp /home/alignment/*.fastq.gz $OUTPUT_DIR && gcloud compute ssh --zone $MACHINE_ZONE $MACHINE_NAME -- gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp ${OUT_BASE}.list $OUTPUT_DIR 
      EXIT_STATUS=$?
      echo "[$(date)] Upload exit status: $EXIT_STATUS"
      echo "[$(date)] Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"
    elif [[ $FAILED_CONTAINER_POLL_COUNT == 5 && $RETRY_COUNTER -lt 4 ]]
    then
      echo "[$(date)] Lost communication with $MACHINE_NAME"
      EXIT_STATUS=-1
    fi
  fi
fi

# Machine may exist even if "docker-machine create" fails.
echo "[$(date)] Cleaning up ..."
gcloud compute instances delete --zone $MACHINE_ZONE --quiet $MACHINE_NAME

exit $EXIT_STATUS
