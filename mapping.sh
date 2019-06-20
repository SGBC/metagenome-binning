#! /usr/bin/env bash

#$ -N ISS reads binning
#$ -M antoine.druart@slu.se
#$ -m seab
#$ -cwd   #Use the directory you're running from
#$ -l h_rt=8:00:0,h_vmem=4G   #Setting running time in hours:min:sec and the memory required for the job
#$ -j y   #Joining the output from standard out and standard error to one file
#$ -pe smp 6-24   #Setting the number of threads for the job to best fit for the system between 1 and N.
#$ -e binning-metrics.log
#$ -o binning-metrics.log
#$ -V

module load bowtie
module load samtools

mkdir -p samples/mapping
bowtie2-build samples/contigs/final.contigs.fa samples/mapping/bt2_index_base
bowtie2 -x samples/mapping/bt2_index_base -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq | samtools view -bS -o samples/mapping/reads_to_sort.bam
samtools sort samples/mapping/reads_to_sort.bam -o samples/mapping/reads.bam
samtools index samples/mapping/reads.bam
