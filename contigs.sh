#! /usr/bin/env bash

#$ -N Megahit-binning
#$ -M antoine.druart@slu.se
#$ -m seab
#$ -cwd   #Use the directory you're running from
#$ -l h_rt=8:00:0,h_vmem=4G   #Setting running time in hours:min:sec and the memory required for the job
#$ -j y   #Joining the output from standard out and standard error to one file
#$ -pe smp 6-24   #Setting the number of threads for the job to best fit for the system between 1 and N.
#$ -e binning-metrics.log
#$ -o binning-metrics.log
#$ -V

module load megahit

rm -rf samples/contigs
megahit -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq -o samples/contigs/
