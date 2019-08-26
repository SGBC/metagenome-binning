#! /usr/bin/env bash

echo "GMM Test with palindrom, any .bam file and .gff file"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall -o results/GMM/F1000_ball_hall
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn -o results/GMM/F1000_dunn
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette -o results/GMM/F1000_silhouette
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index -o results/GMM/F1000_ch-index
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index -o results/GMM/F1000_db-index
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic -o results/GMM/F1000_aic
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic -o results/GMM/F1000_bic
echo "GMM Test without palindrom, any .bam file and .gff file"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal -o results/GMM/F1000_ball_hall_nopal
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal -o results/GMM/F1000_dunn_nopal
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal -o results/GMM/F1000_silhouette_nopal
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal -o results/GMM/F1000_ch-index_nopal
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal -o results/GMM/F1000_db-index_nopal
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic --nopal -o results/GMM/F1000_aic_nopal
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic --nopal -o results/GMM/F1000_bic_nopal

echo "GMM Test with palindrom, any .gff file, with .bam file"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --bam samples/mapping/reads.bam -o results/GMM/F1000_ball_hall_bam
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --bam samples/mapping/reads.bam -o results/GMM/F1000_dunn_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --bam samples/mapping/reads.bam -o results/GMM/F1000_silhouette_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --bam samples/mapping/reads.bam -o results/GMM/F1000_ch-index_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --bam samples/mapping/reads.bam -o results/GMM/F1000_db-index_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic --bam samples/mapping/reads.bam -o results/GMM/F1000_aic_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic --bam samples/mapping/reads.bam -o results/GMM/F1000_bic_bam
echo "GMM Test without palindrom, any .gff file, with .bam file"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --bam samples/mapping/reads.bam -o results/GMM/F1000_ball_hall_nopal_bam
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --bam samples/mapping/reads.bam -o results/GMM/F1000_dunn_nopal_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --bam samples/mapping/reads.bam -o results/GMM/F1000_silhouette_nopal_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --bam samples/mapping/reads.bam  -o results/GMM/F1000_ch-index_nopal_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --bam samples/mapping/reads.bam -o results/GMM/F1000_db-index_nopal_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic --nopal --bam samples/mapping/reads.bam -o results/GMM/F1000_aic_nopal_bam
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic --nopal --bam samples/mapping/reads.bam -o results/GMM/F1000_bic_nopal_bam

echo "GMM Test with palindrom, any .bam file, with .gff file"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ball_hall_cd
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_dunn_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_silhouette_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ch-index_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_db-index_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_aic_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_bic_cd
echo "GMM Test without palindrom, any .bam file, with .gff file"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ball_hall_nopal_cd
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_dunn_nopal_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_silhouette_nopal_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ch-index_nopal_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_db-index_nopal_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic --nopal --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_aic_nopal_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic --nopal --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_bic_nopal_cd

echo "GMM Test with palindrom, with .bam and .gff files"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ball_hall_bam_cd
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_dunn_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_silhouette_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ch-index_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_db-index_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_aic_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_bic_bam_cd
echo "GMM Test without palindrom, with .bam and .gff files"
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ball_hall_nopal_bam_cd
# ./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_dunn_nopal_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_silhouette_nopal_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_ch-index_nopal_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_db-index_nopal_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m aic --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_aic_nopal_bam_cd
./scripts/exploration/gmm.py -i samples/contigs/final.contigs.fa -f 1000 -m bic --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/GMM/F1000_bic_nopal_bam_cd
