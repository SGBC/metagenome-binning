#! /usr/bin/env bash

echo "DBSCAN Test with palindrom, any .bam file and .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall -o results/DBSCAN/F1000_ball_hall
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn -o results/DBSCAN/F1000_dunn
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette -o results/DBSCAN/F1000_silhouette
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index -o results/DBSCAN/F1000_ch-index
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index -o results/DBSCAN/F1000_db-index

echo "DBSCAN Test without palindrom, any .bam file and .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal -o results/DBSCAN/F1000_ball_hall_nopal
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal -o results/DBSCAN/F1000_dunn_nopal
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal -o results/DBSCAN/F1000_silhouette_nopal
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal -o results/DBSCAN/F1000_ch-index_nopal
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal -o results/DBSCAN/F1000_db-index_nopal

echo "DBSCAN Test with palindrom, any .gff file, with .bam file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_ball_hall_bam
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_dunn_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_silhouette_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_ch-index_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_db-index_bam

echo "DBSCAN Test without palindrom, any .gff file, with .bam file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_ball_hall_nopal_bam
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_dunn_nopal_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_silhouette_nopal_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --bam samples/mapping/reads.bam  -o results/DBSCAN/F1000_ch-index_nopal_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1000_db-index_nopal_bam

echo "DBSCAN Test with palindrom, any .bam file, with .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ball_hall_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_dunn_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_silhouette_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ch-index_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_db-index_cd

echo "DBSCAN Test without palindrom, any .bam file, with .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ball_hall_nopal_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_dunn_nopal_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_silhouette_nopal_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ch-index_nopal_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_db-index_nopal_cd


echo "DBSCAN Test with palindrom, with .bam and .gff files"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ball_hall_bam_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_dunn_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_silhouette_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ch-index_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_db-index_bam_cd

echo "DBSCAN Test without palindrom, with .bam and .gff files"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ball_hall_nopal_bam_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_dunn_nopal_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_silhouette_nopal_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ch-index_nopal_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_db-index_nopal_bam_cd



echo "DBSCAN Test with palindrom, any .bam file and .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ball-hall -o results/DBSCAN/F1000_ball_hall_weight
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m dunn -o results/DBSCAN/F1000_dunn_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m silhouette -o results/DBSCAN/F1000_silhouette_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ch-index -o results/DBSCAN/F1000_ch-index_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m db-index -o results/DBSCAN/F1000_db-index_weight

echo "DBSCAN Test without palindrom, any .bam file and .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ball-hall --nopal -o results/DBSCAN/F1000_ball_hall_nopal_weight
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m dunn --nopal -o results/DBSCAN/F1000_dunn_nopal_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m silhouette --nopal -o results/DBSCAN/F1000_silhouette_nopal_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ch-index --nopal -o results/DBSCAN/F1000_ch-index_nopal_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m db-index --nopal -o results/DBSCAN/F1000_db-index_nopal_weight

echo "DBSCAN Test with palindrom, with .bam and .gff files"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ball-hall --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ball_hall_bam_cd_weight
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_dunn_bam_cd_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m silhouette --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_silhouette_bam_cd_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ch-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ch-index_bam_cd_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m db-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_db-index_bam_cd_weight

echo "DBSCAN Test without palindrom, with .bam and .gff files"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ball-hall --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ball_hall_nopal_bam_cd_weight
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m dunn --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_dunn_nopal_bam_cd_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m silhouette --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_silhouette_nopal_bam_cd_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m ch-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_ch-index_nopal_bam_cd_weight
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1000 -w -m db-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/DBSCAN/F1000_db-index_nopal_bam_cd_weight
