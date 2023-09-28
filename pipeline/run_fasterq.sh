#!/bin/bash -e

# This script expects a directory matching the project name
# to be present at /home/blekhman/shared/compendium/accession_lists/
# with a file in it called "SraAccList.txt" listing all the samples
# to be downloaded from SRA.

if [ -z "${PROJECT}" ]
then
      echo "PROJECT variable is EMPTY. Bailing"
      exit 1
fi

echo "${PROJECT},download,start" >> /home/blekhman/shared/compendium/activity.csv

cp -R /home/blekhman/shared/compendium/accession_lists/${PROJECT} /scratch.global/rabdill/bulk/${PROJECT}
cd /scratch.global/rabdill/bulk/${PROJECT}
mkdir -p temp
mkdir -p results
mkdir -p intermediate

echo "Fetching SRA files"
num=1
issues=0
for sample in $(cat SraAccList.txt)
do
    echo "On sample ${num}: $sample"
    problem=0
    if [[ -f "/scratch.global/rabdill/bulk/${PROJECT}/fastq/${sample}_1.fastq" || -f "/scratch.global/rabdill/bulk/${PROJECT}/fastq/${sample}.fastq" ]]; then
        echo "Already found! Continuing."
    else
        /home/blekhman/shared/compendium/code/sratoolkit.2.11.3-centos_linux64/bin/fasterq-dump $sample -e 4 -S -O /scratch.global/rabdill/bulk/${PROJECT}/fastq || problem=1
    fi

    if [[ $problem -gt 0 ]]; then
        echo "SOMETHING WENT WRONG HERE. Continuing..."
        issues=$(($issues + 1))
    fi
    if [[ $issues -gt 10 ]]; then
        echo "Too many issues with this project. Bailing."
        exit 1
    fi
    num=$(($num + 1))
done

cd fastq
ls *_1.fastq | cut -f1 -d "_" > ../samples.txt

# If the samples are only single-ended FASTQ, the format will be different:
if [[ $(cat ../samples.txt | wc -l) -lt 3 ]]; then
    echo "FASTQ files don't end with _1 and _2."
    for f in *.fastq; do
        fnew=`echo $f | sed 's/.fastq/_1.fastq/'`
        mv $f $fnew
    done
fi

# build the sample list again
ls *_1.fastq | cut -f1 -d "_" > ../samples.txt

if [[ $(cat ../samples.txt | wc -l) -lt 3 ]]; then
    echo "Not enough samples downloaded."
    exit 1
fi

echo "${PROJECT},download,end" >> /home/blekhman/shared/compendium/activity.csv

echo "DONE. Submitting dada job."

cd /home/blekhman/shared/compendium/logs

bash ../code/run_dada.slurm ${PROJECT}
