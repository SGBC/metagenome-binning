blk = \e[0m
red = \e[1m\e[31m
green = \e[1m\e[32m
yellow = \e[1m\e[38;5;208m
lblue = \e[1m\e[1;36m

info:
	@echo -e '"make info" \tto show this message \t\t$(green)(implemented)$(blk)'
	@echo -e '"make samples" \tto generate contigs \t\t$(green)(implemented)$(blk)'
	@echo -e "include "make dl_samples", "make extract", "make reads", "make contigs"
	@echo -e '"make bins" \tto generate bins \t\t$(green)(implemented)$(blk)'
	@echo -e "include "make metabat", "make concoct"
	@echo -e '"make metrics"\tto evaluate binning\t\t$(yellow)(not yet fully implemented)$(blk)'
	@echo -e '"make all" \tto generate contigs and bins \t$(yellow)(not yet fully implemented)$(blk)'

samples: dl_samples extract reads contigs

bins: metabat concoct

metrics:
	@echo -e '$(yellow)! Metrics are not yet fully implemented$(blk)'
	

# samples/diamond_db/*: samples/prot_map/*
# 	scripts/diamond_db.py

# samples/prot_map/*: samples/contigs/final.contigs.fa
# 	scripts/prot_map_ref.py

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

reads:
	@echo -e "$(lblue)# Generating reads$(blk)"
	@qsub reads.sh
contigs:
	@echo -e "$(lblue)# Generating contigs$(blk)"
	@qsub contigs.sh

samples/mapping/reads.bam: samples/contigs/final.contigs.fa
	@echo -e "$(lblue)# Mapping reads$(blk)"
	@qsub mapping.sh

metabat: samples/mapping/reads.bam 
	@echo -e "$(yellow)=== Metabat ===$(blk)"
	@qsub metabat.sh

concoct: samples/mapping/reads.bam
	@echo -e "$(yellow)=== Concoct ===$(blk)"
	@qsub concoct.sh

clustering:
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/Optimale10 --bam samples/mapping/reads.bam -c 10
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/ball-hall-30 -m ball-hall --bam samples/mapping/reads.bam -c 30
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/dunn-30 -m dunn --bam samples/mapping/reads.bam -c 30
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/silhouette-30 -m silhouette --bam samples/mapping/reads.bam -c 30
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/ch-index-30 -m ch-index --bam samples/mapping/reads.bam -c 30
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/db-index-30 -m ball-hall --bam samples/mapping/reads.bam -c 30
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/ball-hall-30-nopal -m ball-hall --bam samples/mapping/reads.bam -c 30 --nopal
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/dunn-30-nopal -m dunn --bam samples/mapping/reads.bam -c 30 --nopal
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/silhouette-30-nopal -m silhouette --bam samples/mapping/reads.bam -c 30 --nopal
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/ch-index-30-nopal -m ch-index --bam samples/mapping/reads.bam -c 30 --nopal
	./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -o results/K-means/db-index-30-nopal -m ball-hall --bam samples/mapping/reads.bam -c 30 --nopal
	
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/Optimale10 -c 10
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/ball-hall -m ball-hall
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/dunn -m dunn
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/silhouette -m silhouette
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/ch-index -m ch-index
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/db-index -m db-index
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/ball-hall-nopal -m ball-hall --nopal
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/dunn-nopal -m dunn --nopal
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/silhouette-nopal -m silhouette --nopal
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/ch-index-nopal -m ch-index --nopal
	./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --bam samples/mapping/reads.bam -o results/H_clust/db-index-nopal -m db-index --nopal

report:
	scripts/validation/fast_metrics.py -i results/metabat/fasta_bins/* -m samples/seq_index.json -n "Metabat" -o graphs/metabat/ --map
	scripts/validation/fast_metrics.py -i results/concoct/fasta_bins/* -m samples/seq_index.json -n "concoct" -o graphs/concoct/ --map

	scripts/validation/fast_metrics.py -i results/H_clust/Optimale10/* -m samples/seq_index.json -n "H_clust Optimal (nb_clust = 10)" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/ball-hall/* -m samples/seq_index.json -n "H_clust ball-hall" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/ball-hall-nopal/* -m samples/seq_index.json -n "H_clust ball-hall nopal" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/ch-index/* -m samples/seq_index.json -n "H_clust ch-index" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/ch-index-nopal/* -m samples/seq_index.json -n "H_clust ch-index-nopal" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/db-index/* -m samples/seq_index.json -n "H_clust db-index" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/db-index-nopal/* -m samples/seq_index.json -n "H_clust db-index-nopal" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/dunn/* -m samples/seq_index.json -n "H_clust dunn" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/dunn-nopal/* -m samples/seq_index.json -n "H_clust dunn-nopal" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/silhouette/* -m samples/seq_index.json -n "H_clust silhouette" -o graphs/H_clust/ --map
	scripts/validation/fast_metrics.py -i results/H_clust/silhouette-nopal/* -m samples/seq_index.json -n "H_clust silhouette-nopal" -o graphs/H_clust/ --map

	scripts/validation/fast_metrics.py -i results/K-means/Optimale10/* -m samples/seq_index.json -n "K-means Optimal (nb_clust = 10)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/ball-hall-30/* -m samples/seq_index.json -n "K-means ball-hall (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/ball-hall-30-nopal/* -m samples/seq_index.json -n "K-means ball-hall nopal (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/ch-index-30/* -m samples/seq_index.json -n "K-means ch-index (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/ch-index-30-nopal/* -m samples/seq_index.json -n "K-means ch-index-nopal (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/db-index-30/* -m samples/seq_index.json -n "K-means db-index (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/db-index-30-nopal/* -m samples/seq_index.json -n "K-means db-index-nopal (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/dunn-30/* -m samples/seq_index.json -n "K-means dunn (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/dunn-30-nopal/* -m samples/seq_index.json -n "K-means dunn-nopal (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/silhouette-30/* -m samples/seq_index.json -n "K-means silhouette (nb_max_clust = 29)" -o graphs/K-means/ --map
	scripts/validation/fast_metrics.py -i results/K-means/silhouette-30-nopal/* -m samples/seq_index.json -n "K-means silhouette-nopal (nb_max_clust = 29)" -o graphs/K-means/ --map
