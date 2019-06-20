#! /usr/bin/env bash

#$ -N metabat-binning
#$ -M antoine.druart@slu.se
#$ -m seab
#$ -cwd   #Use the directory you're running from
#$ -l h_rt=8:00:0,h_vmem=4G   #Setting running time in hours:min:sec and the memory required for the job
#$ -j y   #Joining the output from standard out and standard error to one file
#$ -pe smp 6-24   #Setting the number of threads for the job to best fit for the system between 1 and N.
#$ -e binning-metrics.log
#$ -o binning-metrics.log
#$ -V

module load metabat
module load python/3.7.3

runMetaBat.sh -m 1500 samples/contigs/final.contigs.fa samples/mapping/reads.bam
mv final.contigs.fa.metabat-bins1500 samples/metabat/
mv final.contigs.fa.depth.txt samples/metabat/
mv final.contigs.fa.paired.txt samples/metabat/
mv samples/metabat/final.contigs.fa.metabat-bins1500 samples/metabat/fasta_bins
scripts/metabin_rename.py