blk = \e[0m
red = \e[1m\e[31m
green = \e[1m\e[32m
yellow = \e[1m\e[38;5;208m
lblue = \e[38;5;50m

info:
	@echo -e '"make info" \tto show this message \t\t$(green)(implemented)$(blk)'
	@echo -e '"make samples" \tto generate contigs \t\t$(green)(implemented)$(blk)'
	@echo -e '"make bins" \tto generate bins \t\t$(yellow)(partially implemented)$(blk)'
	@echo -e '"make all" \tto generate contigs and bins \t$(yellow)(partially implemented)$(blk)'

samples: dl_samples extract sort_entry reads contigs

bins: metabat
	@echo -e '$(yellow)! The feature $(blk)$(lblue)"bins"$(yellow) is partially implemented and in testing :s$(blk)'

all:
	samples bins

dl_samples:
	@mkdir -p samples
	@echo "$(lblue)# Download samples(10)$(blk)"
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
	@mkdir -p samples/archives
	@mkdir -p samples/complete_genomes
	@echo "$(lblue)# Extract the data$(blk)"
	@gunzip -kN samples/*.gz
	@echo "$(lblue)# Organize the files$(blk)"
	@mv samples/*.gz samples/archives/
	@mv samples/*.fna samples/complete_genomes
	@echo It looks clean :D

sort_entry:
	@echo "$(lblue)# Running chroplasmitor (sorts entries)$(blk)"
	@scripts/chroplasmitor.py -g samples/complete_genomes/*.fna -o samples

reads:
	@echo "$(lblue)# Generating reads$(blk)"
	@mkdir -p samples/reads
	@iss generate --genomes samples/chromosomes/*.fna --abundance_file scripts/abundance_file --model hiseq --output samples/reads/reads

contigs:
	@echo "$(lblue)# Generating contigs$(blk)"
	@megahit -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq -o samples/contigs/

mapping:
	@echo "$(lblue)# Mapping reads$(blk)"
	@mkdir -p samples/metabat
	@mkdir -p samples/metabat/map
	@echo "$(lblue)# Run bowtie2-build$(blk)"
	@bowtie2-build samples/contigs/final.contigs.fa samples/metabat/map/bt2_index_base
	@echo "$(lblue)# Run bowtie2 and samtools view$(blk)"
	@bowtie2 -x samples/metabat/map/bt2_index_base -1 samples/reads/reads_R1.fastq -2 samples/reads/reads_R2.fastq | samtools view -bS -o samples/metabat/map/reads_to_sort.bam
	@echo "$(lblue)# Run samtools sort$(blk)"
	@samtools sort samples/metabat/map/reads_to_sort.bam -o reads.bam
	@echo "$(lblue)# Run samtools index$(blk)"
	@samtools index reads.bam
	@mv reads.bam samples/metabat/
	@mv reads.bam.bai samples/metabat/

metabat: mapping
	@echo "$(lblue)# Run metabat$(blk)"
	@runMetaBat.sh -m 1500 samples/contigs/final.contigs.fa samples/metabat/reads.bam
	@rm -rf samples/metabat/final.contigs.fa.metabat-bins1500
	@mv final.contigs.fa.metabat-bins1500 samples/metabat/
	@mv final.contigs.fa.depth.txt samples/metabat/
	@mv final.contigs.fa.paired.txt samples/metabat/

concoct:
	@echo "$(yellow)! Concoct implementation in progress$(blk)"
	@mkdir -p samples/concoct
	@mkdir -p samples/concoct/map
	@source activate concoct_env
	@bowtie2-build samples/contigs/final.contigs.fa samples/concoct/map/bt2_index_base
