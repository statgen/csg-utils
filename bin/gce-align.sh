#!/bin/bash
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
  1) BWA_PARIED_OPT="-p" ;;
esac

MACHINE_NAME="align-${SAMPLE_ID}-${RG_IDENT}"
MACHINE_TAG=$(basename $1 .fastq.gz | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)
MACHINE_TYPE_OPTS="--custom-cpu 32 --custom-memory 64GiB"
#MACHINE_TYPE_OPTS="--machine-type  n1-standard-32"

gcloud compute instances create --tags $MACHINE_TAG --image ubuntu-1404-alignment $MACHINE_TYPE_OPTS --boot-disk-size 300 --preemptible $MACHINE_NAME

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
      echo "Uploading "$f" ..."
      START_TIME=$(date +%s)
      gcloud compute copy-files $f $MACHINE_NAME":/home/alignment/input.fastq.gz"
      EXIT_STATUS=$?
      echo "Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"

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
            
            echo "Fetching logs ..."
            gcloud compute ssh $MACHINE_NAME --command "sudo docker logs $CONTAINER_ID"
            if [[ $? != 0 ]]
            then
              echo "Fetching logs FAILED!"
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
            echo "Failed Count: "$FAILED_CONTAINER_POLL_COUNT
          fi
        done

        if [[ $CONTAINER_IS_RUNNING == 0 && $EXIT_STATUS == 0 ]]
        then
          OUTPUT_FILE=$OUTPUT_DIR"/"$(basename $f .fastq.gz)".cram"
          echo "Downloading "$OUTPUT_FILE" ..."
          START_TIME=$(date +%s)
          gcloud compute copy-files $MACHINE_NAME":/home/alignment/output.cram" $OUTPUT_FILE && gcloud compute copy-files $MACHINE_NAME":/home/alignment/output.cram.ok" $OUTPUT_FILE".ok"
          EXIT_STATUS=$?
          echo 'Download exit status: '$EXIT_STATUS
          echo "Elapsed time: "$(( $(date +%s) - $START_TIME ))"s"
          break
        elif [[ $FAILED_CONTAINER_POLL_COUNT == 5 && $RETRY_COUNTER -lt 4 ]]
        then
          #sleep 300s
          echo "Machine stopped. Restarting "$MACHINE_NAME" ..."
          gcloud compute instances start $MACHINE_NAME
        fi

        let RETRY_COUNTER++
      done
    done
  fi
fi

# Machine may exist even if "docker-machine create" fails.
echo 'Cleaning up ...'
gcloud compute instances delete --quiet $MACHINE_NAME

exit $EXIT_STATUS
