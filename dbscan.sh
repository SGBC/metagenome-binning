#! /usr/bin/env bash

echo "DBSCAN Test with palindrom, any .bam file and .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall -o results/DBSCAN/F1500_ball_hall
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn -o results/DBSCAN/F1500_dunn
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette -o results/DBSCAN/F1500_silhouette
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index -o results/DBSCAN/F1500_ch-index
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index -o results/DBSCAN/F1500_db-index

echo "DBSCAN Test without palindrom, any .bam file and .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall --nopal -o results/DBSCAN/F1500_ball_hall_nopal
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn --nopal -o results/DBSCAN/F1500_dunn_nopal
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette --nopal -o results/DBSCAN/F1500_silhouette_nopal
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index --nopal -o results/DBSCAN/F1500_ch-index_nopal
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index --nopal -o results/DBSCAN/F1500_db-index_nopal

echo "DBSCAN Test with palindrom, any .gff file, with .bam file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_ball_hall_bam
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_dunn_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_silhouette_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_ch-index_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_db-index_bam

echo "DBSCAN Test without palindrom, any .gff file, with .bam file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_ball_hall_nopal_bam
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_dunn_nopal_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_silhouette_nopal_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index --nopal --bam samples/mapping/reads.bam  -o results/DBSCAN/F1500_ch-index_nopal_bam
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index --nopal --bam samples/mapping/reads.bam -o results/DBSCAN/F1500_db-index_nopal_bam

echo "DBSCAN Test with palindrom, any .bam file, with .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ball_hall_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_dunn_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_silhouette_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ch-index_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_db-index_cd

echo "DBSCAN Test without palindrom, any .bam file, with .gff file"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall --nopal --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ball_hall_nopal_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn --nopal --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_dunn_nopal_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette --nopal --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_silhouette_nopal_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index --nopal --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ch-index_nopal_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index --nopal --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_db-index_nopal_cd


echo "DBSCAN Test with palindrom, with .bam and .gff files"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ball_hall_bam_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_dunn_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_silhouette_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ch-index_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_db-index_bam_cd

echo "DBSCAN Test without palindrom, with .bam and .gff files"
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ball-hall --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ball_hall_nopal_bam_cd
# ./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m dunn --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_dunn_nopal_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m silhouette --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_silhouette_nopal_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m ch-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_ch-index_nopal_bam_cd
./scripts/exploration/dbscan.py -i samples/contigs/final.contigs.fa -f 1500 -m db-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_prot.gff -o results/DBSCAN/F1500_db-index_nopal_bam_cd
