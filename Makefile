blk = \e[0m
red = \e[1m\e[31m
green = \e[1m\e[32m
yellow = \e[1m\e[38;5;208m
lblue = \e[1m\e[1;36m

info:
	@echo -e '"make info" \tto show this message \t\t$(green)(implemented)$(blk)'
	@echo -e '"make samples" \tto generate contigs \t\t$(green)(implemented)$(blk)'
	@echo -e 'include "make dl_samples", "make extract", "make reads", "make contigs"'
	@echo -e '"make bins" \tto generate bins \t\t$(green)(implemented)$(blk)'
	@echo -e 'include "make metabat", "make concoct"'
	@echo -e '"make build_metrics_index"\tto compute an alignement between contigs and originals chomosomes.\t$(green)(implemented)$(blk)'
	@echo -e '"make all" \tto generate contigs, metrics_index and bins by concoct, metabat and experimentals algorithme \t$(yellow)(not yet fully implemented)$(blk)'

samples: dl_samples extract reads contigs

bins: metabat concoct

build_metrics_index: samples/diamond_db/*
	#@echo -e '$(yellow)Metrics are not yet fully implemented$(blk)'
	@scripts/validation/precomputed_refindex.py -i samples/contigs/final.contigs.fa -o samples/ -r samples/chromosomes/all_chromo.fna --db samples/diamond_db/* --index samples/chromosomes/index_chromo

samples/diamond_db/*: samples/prot_map/[A-Z]*.faa
	scripts/tools/diamond_db.py -i samples/prot_map/*.faa 

samples/prot_map/*: samples/contigs/final.contigs.fa
	scripts/tools/prot_map_ref.py -i samples/chromosomes/[A-Z]*.fna

all: samples build_metrics_index bins 

dl_samples:
	@mkdir -p samples
	@echo -e "$(lblue)# Download samples(10)$(blk)"
	@echo Download Bacillus subtilis
	@wget -qO samples/Bacillus_subtilis.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz
	@echo Download Cryptococcus neoformans
	@wget -qO samples/Cryptococcus_neoformans.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz
	@echo Download Enterococcus faecalis
	@wget -qO samples/Enterococcus_faecalis.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/785/GCF_000007785.1_ASM778v1/GCF_000007785.1_ASM778v1_genomic.fna.gz
	@echo Download Escherichia coli
	@wget -qO samples/Escherichia_coli.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
	@echo Download Lactobacillus fermentum
	@wget -qO samples/Lactobacillus_fermentum.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/145/GCF_000010145.1_ASM1014v1/GCF_000010145.1_ASM1014v1_genomic.fna.gz
	@echo Download Listeria monocytogenes
	@wget -qO samples/Listeria_monocytogenes.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/035/GCF_000196035.1_ASM19603v1/GCF_000196035.1_ASM19603v1_genomic.fna.gz
	@echo Download Pseudomonas aeruginosa
	@wget -qO samples/Pseudomonas_aeruginosa.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz
	@echo Download Saccharomyces cerevisiae
	@wget -qO samples/Saccharomyces_cerevisiae.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
	@echo Download Salmonella enterica
	@wget -qO samples/Salmonella_enterica.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/995/GCF_000195995.1_ASM19599v1/GCF_000195995.1_ASM19599v1_genomic.fna.gz
	@echo Download Staphylococcus aureus
	@wget -qO samples/Staphylococcus_aureus.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz

extract:
	@mkdir -p samples/complete_genomes
	@echo -e "$(lblue)# Extract the data$(blk)"
	@gunzip -Nf samples/*.gz
	@echo -e "$(lblue)# Organize the files$(blk)"
	@mv samples/*.fna samples/complete_genomes
	@echo It looks clean :D

reads: 
	@echo -e "$(lblue)# Generating reads$(blk)"
	@./scripts/tools/chroplasmitor.py -g samples/complete_genomes/*.fna -o samples
	@./scripts/tools/chromo_compiler.py
	@mkdir -p samples/reads
	@iss generate --genomes samples/chromosomes/all_chromo.fna --abundance_file scripts/abundance_file --model hiseq --output samples/reads/reads --n_reads 17M --cpus `grep -c ^processor /proc/cpuinfo`

samples/reads/reads_R1.fastq: samples/complete_genomes/Bacillus_subtilis.fna samples/complete_genomes/Cryptococcus_neoformans.fna samples/complete_genomes/Enterococcus_faecalis.fna samples/complete_genomes/Escherichia_coli.fna samples/complete_genomes/Lactobacillus_fermentum.fna samples/complete_genomes/Listeria_monocytogenes.fna samples/complete_genomes/Pseudomonas_aeruginosa.fna samples/complete_genomes/Saccharomyces_cerevisiae.fna samples/complete_genomes/Salmonella_enterica.fna samples/complete_genomes/Staphylococcus_aureus.fna
	make reads

samples/reads/reads_R2.fastq: samples/complete_genomes/Bacillus_subtilis.fna samples/complete_genomes/Cryptococcus_neoformans.fna samples/complete_genomes/Enterococcus_faecalis.fna samples/complete_genomes/Escherichia_coli.fna samples/complete_genomes/Lactobacillus_fermentum.fna samples/complete_genomes/Listeria_monocytogenes.fna samples/complete_genomes/Pseudomonas_aeruginosa.fna samples/complete_genomes/Saccharomyces_cerevisiae.fna samples/complete_genomes/Salmonella_enterica.fna samples/complete_genomes/Staphylococcus_aureus.fna
	make reads

contigs: 
	@echo -e "$(lblue)# Generating contigs$(blk)"
	@rm -rf samples/contigs
	@megahit -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq -o samples/contigs/

samples/mapping/reads.bam: samples/contigs/final.contigs.fa
	@echo -e "$(lblue)# Mapping reads$(blk)"
	@mkdir -p samples/mapping
	echo bowtie2-build
	@bowtie2-build samples/contigs/final.contigs.fa samples/mapping/bt2_index_base
	echo bowtie2 -x
	@bowtie2 -p `grep -c ^processor /proc/cpuinfo` -x samples/mapping/bt2_index_base -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq | samtools view -bS -o samples/mapping/reads_to_sort.bam
	echo samtools
	@samtools sort samples/mapping/reads_to_sort.bam -o samples/mapping/reads.bam
	@samtools index samples/mapping/reads.bam

samples/prot_map/contigs_genes.gff: samples/contigs/final.contigs.fa
	mkdir -p samples/prot_map
	prodigal -i samples/contigs/final.contigs.fa -g 11 -f gff -m -o samples/prot_map/contigs_genes.gff

samples/contigs/final.contigs.fa: samples/reads/reads_R1.fastq samples/reads/reads_R2.fastq
	make contigs

results/metabat/fasta_bins/*: samples/mapping/reads.bam
	@echo -e "$(yellow)=== Metabat ===$(blk)"
	@runMetaBat.sh -m 1500 samples/contigs/final.contigs.fa samples/mapping/reads.bam
	@mkdir -p results/metabat/fasta_bins
	@mv final.contigs.fa.metabat-bins1500 results/metabat/
	@mv final.contigs.fa.depth.txt results/metabat/
	@mv final.contigs.fa.paired.txt results/metabat/
	@mv results/metabat/final.contigs.fa.metabat-bins1500/* results/metabat/fasta_bins
	@rm -r results/metabat/final.contigs.fa.metabat-bins1500
	@scripts/tools/metabin_rename.py

results/concoct/fasta_bins/*: samples/mapping/reads.bam
	@echo -e "$(yellow)=== Concoct ===$(blk)"
	@mkdir -p results/concoct
	@echo cut_up_fasta
	@cut_up_fasta.py samples/contigs/final.contigs.fa  -c 10000 -o 0 --merge_last -b results/concoct/contigs_10K.bed > results/concoct/contigs_10K.fa
	@echo concoct_coverage_table
	@concoct_coverage_table.py results/concoct/contigs_10K.bed samples/mapping/reads.bam > results/concoct/coverage_table.tsv
	@concoct --composition_file results/concoct/contigs_10K.fa --coverage_file results/concoct/coverage_table.tsv -b results/concoct/output/
	@rm -rf results/concoct/fasta_bins
	@merge_cutup_clustering.py results/concoct/output/clustering_gt1000.csv > results/concoct/output/clustering_merged.csv
	@mkdir results/concoct/fasta_bins
	@extract_fasta_bins.py samples/contigs/final.contigs.fa results/concoct/output/clustering_merged.csv --output_path results/concoct/fasta_bins

metabat: results/metabat/fasta_bins/* 
	@echo -e "$(green)=== Metabat done ===$(blk)"

concoct: results/concoct/fasta_bins/*
	@echo -e "$(green)=== Concoct done ===$(blk)"

clustering: samples/mapping/reads.bam samples/prot_map/contigs_genes.gff
	scripts/exploration/hierarchical_clustering.sh
	scripts/exploration/kmeans.sh
	scripts/exploration/dbscan.sh
	scripts/exploration/gmm.sh
	scripts/exploration/vbgmm.sh
	scripts/exploration/affinity_propagation.sh

report:
	./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/METABAT/ -i results/metabat/fasta_bins/* -n "metabat_set"
	./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/CONCOCT/ -i results/concoct/fasta_bins/* -n "concoct_set"
	scripts/validation/hierarchical_clustering_metrics.sh
	scripts/validation/kmeans_metrics.sh
	scripts/validation/dbscan_metrics.sh
	scripts/validation/gmm_metrics.sh
	scripts/validation/vbgmm_metrics.sh
	scripts/validation/affinity_propagation_metrics.sh
