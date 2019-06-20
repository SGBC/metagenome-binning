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

module load biopython
module load insilicoseq
module load python/3.7.3

scripts/chroplasmitor.py -g samples/complete_genomes/*.fna -o samples
scripts/chromo_compiler.py
mkdir -p samples/reads
iss generate --genomes samples/chromosomes/all_chromo.fna --abundance_file scripts/abundance_file --model hiseq --output samples/reads/reads --n_reads 20M --cpus `grep -c ^processor /proc/cpuinfo`
