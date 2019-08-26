#! /usr/bin/env bash

# echo "K-MEANS Test with palindrom, any .bam file and .gff file Optimum clusters"
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 -o results/K-MEANS/SIMPLE_OPTIMUM
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 --nopal -o results/K-MEANS/SIMPLE_OPTIMUM_nopal
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 --bam samples/mapping/reads.bam -o results/K-MEANS/SIMPLE_OPTIMUM_bam
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 --nopal --bam samples/mapping/reads.bam -o results/K-MEANS/SIMPLE_OPTIMUM_nopal_bam
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/SIMPLE_OPTIMUM_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 --nopal --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/SIMPLE_OPTIMUM_nopal_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/SIMPLE_OPTIMUM_bam_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa  -c 10 --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/SIMPLE_OPTIMUM_nopal_bam_cd

# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 -o results/K-MEANS/F2000_OPTIMUM
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 --nopal -o results/K-MEANS/F2000_OPTIMUM_nopal
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_OPTIMUM_bam
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 --nopal --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_OPTIMUM_nopal_bam
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_OPTIMUM_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 --nopal --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_OPTIMUM_nopal_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_OPTIMUM_bam_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -c 10 --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_OPTIMUM_nopal_bam_cd


echo "K-MEANS Test with palindrom, any .bam file and .gff file"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall -o results/K-MEANS/F2000_ball_hall
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn -o results/K-MEANS/F2000_dunn
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette -o results/K-MEANS/F2000_silhouette
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index -o results/K-MEANS/F2000_ch-index
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index -o results/K-MEANS/F2000_db-index

echo "K-MEANS Test without palindrom, any .bam file and .gff file"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal -o results/K-MEANS/F2000_ball_hall_nopal
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal -o results/K-MEANS/F2000_dunn_nopal
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal -o results/K-MEANS/F2000_silhouette_nopal
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal -o results/K-MEANS/F2000_ch-index_nopal
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal -o results/K-MEANS/F2000_db-index_nopal

echo "K-MEANS Test with palindrom, any .gff file, with .bam file"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_ball_hall_bam
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_dunn_bam
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_silhouette_bam
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_ch-index_bam
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_db-index_bam

echo "K-MEANS Test without palindrom, any .gff file, with .bam file"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_ball_hall_nopal_bam
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_dunn_nopal_bam
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_silhouette_nopal_bam
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal --bam samples/mapping/reads.bam  -o results/K-MEANS/F2000_ch-index_nopal_bam
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal --bam samples/mapping/reads.bam -o results/K-MEANS/F2000_db-index_nopal_bam

echo "K-MEANS Test with palindrom, any .bam file, with .gff file"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ball_hall_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_dunn_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_silhouette_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ch-index_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_db-index_cd

echo "K-MEANS Test without palindrom, any .bam file, with .gff file"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ball_hall_nopal_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_dunn_nopal_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_silhouette_nopal_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ch-index_nopal_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_db-index_nopal_cd


echo "K-MEANS Test with palindrom, with .bam and .gff files"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ball_hall_bam_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_dunn_bam_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_silhouette_bam_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ch-index_bam_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_db-index_bam_cd

echo "K-MEANS Test without palindrom, with .bam and .gff files"
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ball_hall_nopal_bam_cd
# ./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_dunn_nopal_bam_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_silhouette_nopal_bam_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_ch-index_nopal_bam_cd
./scripts/exploration/k-means.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/K-MEANS/F2000_db-index_nopal_bam_cd
