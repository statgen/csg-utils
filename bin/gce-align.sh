#!/bin/bash
RG_LINE=$1
shift
OUTPUT_DIR=$1
shift
INPUT_FILES=$@

set -o pipefail

MACHINE_NAME="align-"$(head -c16 /dev/urandom | xxd -p | tr -d \\n)
MACHINE_TAG=$(basename $1 .fastq.gz | tr "[:upper:]" "[:lower:]" | sed "s/[^a-z0-9]/-/g" | head -c62)
MACHINE_TYPE_OPTS="--custom-cpu 32 --custom-memory 64GiB"
#MACHINE_TYPE_OPTS="--machine-type  n1-standard-32"

gcloud compute instances create --tags $MACHINE_TAG --image ubuntu-1404-alignment $MACHINE_TYPE_OPTS --boot-disk-size 300 --preemptible $MACHINE_NAME

EXIT_STATUS=$?
if [[ $EXIT_STATUS == 0 ]]
then
  #COMMAND="set -o pipefail; bwa mem -t 32 -K 100000000 -Y -p -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samtools sort -@ 32 -m 2000M --reference /home/alignment/ref/hs38DH.fa -O cram -o /home/alignment/output.cram -T /home/alignment/sort.temp -"
  COMMAND="set -o pipefail; rm -f /home/alignment/output.cram /home/alignment/output.cram.ok; bwa mem -t 32 -K 100000000 -Y -p -R '$RG_LINE' /home/alignment/ref/hs38DH.fa /home/alignment/input.fastq.gz | samblaster -a --addMateTags | samtools view -@ 32 -T /home/alignment/ref/hs38DH.fa -C -o /home/alignment/output.cram - && touch /home/alignment/output.cram.ok"
  
  CONTAINER_ID=$(gcloud compute ssh $MACHINE_NAME -- sudo docker create -v "/home/alignment:/home/alignment" statgen/alignment /bin/bash -c \""$COMMAND"\")

  EXIT_STATUS=$?
  if [[ $EXIT_STATUS == 0 ]]
  then  
    for f in $INPUT_FILES
    do
      echo "Uploading "$f" ..."
      gcloud compute copy-files $f $MACHINE_NAME":/home/alignment/input.fastq.gz"
      EXIT_STATUS=$?

      COUNTER=0
      while [[ $EXIT_STATUS == 0 && $COUNTER -lt 5 ]]
      do  

        gcloud compute ssh $MACHINE_NAME -- sudo docker start $CONTAINER_ID

        CONTAINER_IS_RUNNING=1
        MACHINE_IS_RUNNING=1
        while [[ $CONTAINER_IS_RUNNING != 0 && $MACHINE_IS_RUNNING != 0 ]]
        do
          sleep 60s

          if [[ $(gcloud compute ssh $MACHINE_NAME -- sudo docker ps -al --format {{.Status}}) =~ Exited\ \((.*)\) ]]
          then
            CONTAINER_IS_RUNNING=0
            EXIT_STATUS=${BASH_REMATCH[1]}
          fi

          if [[ $(gcloud compute instances list $MACHINE_NAME | grep RUNNING | wc -l) == 0 ]]
          then
            MACHINE_IS_RUNNING=0
          fi

        done

        if [[ $CONTAINER_IS_RUNNING == 0 && $EXIT_STATUS == 0 ]]
        then
          echo "Fetching logs ..."
          gcloud compute ssh $MACHINE_NAME --command "sudo docker logs $CONTAINER_ID"
          if [[ $? != 0 ]]
          then
            echo "Fetching logs FAILED!"
          fi

          OUTPUT_FILE=$OUTPUT_DIR"/"$(basename $f .fastq.gz)".cram"
          echo "Downloading "$OUTPUT_FILE" ..."
          gcloud compute copy-files $MACHINE_NAME":/home/alignment/output.cram" $OUTPUT_FILE && gcloud compute copy-files $MACHINE_NAME":/home/alignment/output.cram.ok" $OUTPUT_FILE".ok"
          EXIT_STATUS=$?
          break
        elif [[ $MACHINE_IS_RUNNING == 0 ]] 
        then
          sleep 300s
          gcloud compute instances start $MACHINE_NAME
        fi

        let COUNTER++
      done
    done
  fi
fi

# Machine may exist even if "docker-machine create" fails.
gcloud compute instances delete --quiet $MACHINE_NAME

exit $EXIT_STATUS
