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

module load concoct


mkdir -p samples/concoct
cut_up_fasta.py samples/contigs/final.contigs.fa  -c 10000 -o 0 --merge_last -b samples/concoct/contigs_10K.bed > samples/concoct/contigs_10K.fa
concoct_coverage_table.py samples/concoct/contigs_10K.bed samples/mapping/reads.bam > samples/concoct/coverage_table.tsv
concoct --composition_file samples/concoct/contigs_10K.fa --coverage_file samples/concoct/coverage_table.tsv -b samples/concoct/output/
rm -rf samples/concoct/fasta_bins
merge_cutup_clustering.py samples/concoct/output/clustering_gt1000.csv > samples/concoct/output/clustering_merged.csv
mkdir samples/concoct/fasta_bins
extract_fasta_bins.py samples/contigs/final.contigs.fa samples/concoct/output/clustering_merged.csv --output_path samples/concoct/fasta_bins