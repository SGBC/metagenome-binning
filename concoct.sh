#! /usr/bin/env bash

#$ -N concoct-binning
#$ -M antoine.druart@slu.se
#$ -m seab
#$ -cwd   #Use the directory you're running from
#$ -l h_rt=8:00:0,h_vmem=4G   #Setting running time in hours:min:sec and the memory required for the job
#$ -j y   #Joining the output from standard out and standard error to one file
#$ -pe smp 6-24   #Setting the number of threads for the job to best fit for the system between 1 and N.
#$ -e binning-metrics.log
#$ -o binning-metrics.log
#$ -V

module load concoct


mkdirs -p results/concoct
cut_up_fasta.py samples/contigs/final.contigs.fa  -c 10000 -o 0 --merge_last -b results/concoct/contigs_10K.bed > results/concoct/contigs_10K.fa
concoct_coverage_table.py results/concoct/contigs_10K.bed samples/mapping/reads.bam > results/concoct/coverage_table.tsv
concoct --composition_file results/concoct/contigs_10K.fa --coverage_file results/concoct/coverage_table.tsv -b results/concoct/output/
rm -rf results/concoct/fasta_bins
merge_cutup_clustering.py results/concoct/output/clustering_gt1000.csv > results/concoct/output/clustering_merged.csv
mkdir results/concoct/fasta_bins
extract_fasta_bins.py results/contigs/final.contigs.fa results/concoct/output/clustering_merged.csv --output_path results/concoct/fasta_bins