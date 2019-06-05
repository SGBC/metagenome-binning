blk = \e[0m
red = \e[1m\e[31m
green = \e[1m\e[32m
yellow = \e[1m\e[38;5;208m
lblue = \e[1m\e[1;36m

info:
	@echo -e '"make info" \tto show this message \t\t$(green)(implemented)$(blk)'
	@echo -e '"make samples" \tto generate contigs \t\t$(green)(implemented)$(blk)'
	@echo -e '"make bins" \tto generate bins \t\t$(green)(implemented)$(blk)'
	@echo -e '"make metrics"\tto evaluate binning\t\t$(yellow)(not yet fully implemented)$(blk)'
	@echo -e '"make server"\tto import all necessary to run the makefile on the HGEN server'
	@echo -e '"make all" \tto generate contigs and bins \t$(yellow)(not yet fully implemented)$(blk)'

samples: dl_samples extract sort_entry reads contigs

server:
	bash scripts/modules_load

bins: metabat #concoct

metrics:
	@echo -e '$(yellow)! Metrics are not yet fully implemented$(blk)'
	@cat samples/metabat/time
	@cat samples/concoct/time


all: samples bins 

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

sort_entry:
	@echo -e "$(lblue)# Running chroplasmitor (sorts entries)$(blk)"
	@scripts/chroplasmitor.py -g samples/complete_genomes/*.fna -o samples

reads:
	@echo -e "$(lblue)# Generating reads$(blk)"
	@mkdir -p samples/reads
	@iss generate --genomes samples/chromosomes/*.fna --abundance_file scripts/abundance_file --model hiseq --output samples/reads/reads --n_reads 8M --cpus `grep -c ^processor /proc/cpuinfo`

contigs:
	@echo -e "$(lblue)# Generating contigs$(blk)"
	@megahit -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq -o samples/contigs/

samples/mapping/reads.bam: samples/contigs/final.contigs.fa
	@echo -e "$(lblue)# Mapping reads$(blk)"
	@mkdir -p samples/mapping
	@echo -e "$(lblue)# Run bowtie2-build$(blk)"
	@bowtie2-build samples/contigs/final.contigs.fa samples/mapping/bt2_index_base
	@echo -e "$(lblue)# Run bowtie2 and samtools view$(blk)"
	@bowtie2 -x samples/mapping/bt2_index_base -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq | samtools view -bS -o samples/mapping/reads_to_sort.bam
	@echo -e "$(lblue)# Run samtools sort$(blk)"
	@samtools sort samples/mapping/reads_to_sort.bam -o samples/mapping/reads.bam
	@echo -e "$(lblue)# Run samtools index$(blk)"
	@samtools index samples/mapping/reads.bam

metabat: samples/mapping/reads.bam 
	@echo -e "$(yellow)=== Metabat ===$(blk)"
	@mkdir -p samples/metabat
	@echo "=== Time running Metabat ===" > samples/metabat/time
	@echo Start @\t `date` >> samples/metabat/time
	@echo -e "$(lblue)# Run metabat$(blk)"
	@runMetaBat.sh -m 1500 samples/contigs/final.contigs.fa samples/mapping/reads.bam
	@rm -rf samples/metabat/final.contigs.fa.metabat-bins1500
	@mv final.contigs.fa.metabat-bins1500 samples/metabat/
	@mv final.contigs.fa.depth.txt samples/metabat/
	@mv final.contigs.fa.paired.txt samples/metabat/
	@mv samples/metabat/final.contigs.fa.metabat-bins1500 samples/metabat/fasta_bins
	@scripts/metabin_rename.py
	@echo End @\t `date` >> samples/metabat/time

concoct: samples/mapping/reads.bam
	@echo -e "$(yellow)=== Concoct ===$(blk)"
	@mkdir -p samples/concoct
	@echo -e "=== Time running Concoct ===" > samples/concoct/time
	@echo Start @\t `date` >> samples/concoct/time
	@echo -e "$(lblue)# Cut contigs in small part$(blk)"
	@cut_up_fasta.py samples/contigs/final.contigs.fa  -c 10000 -o 0 --merge_last -b samples/concoct/contigs_10K.bed > samples/concoct/contigs_10K.fa
	@echo -e "$(lblue)# Generate table of coverage depth$(blk)"
	@concoct_coverage_table.py samples/concoct/contigs_10K.bed samples/mapping/reads.bam > samples/concoct/coverage_table.tsv
	@echo -e "$(lblue)# Run concoct$(blk)"
	@concoct --composition_file samples/concoct/contigs_10K.fa --coverage_file samples/concoct/coverage_table.tsv -b samples/concoct/output/
	@echo -e "$(lblue)# Merge subcontigs clustering$(blk)"
	@rm -rf samples/concoct/fasta_bins
	@merge_cutup_clustering.py samples/concoct/output/clustering_gt1000.csv > samples/concoct/output/clustering_merged.csv
	@mkdir samples/concoct/fasta_bins
	@echo -e "$(lblue)# Extract Bins$(blk)"
	@extract_fasta_bins.py samples/contigs/final.contigs.fa samples/concoct/output/clustering_merged.csv --output_path samples/concoct/fasta_bins
	@echo End @\t `date` >> samples/concoct/time
