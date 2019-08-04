#! /usr/bin/env bash

echo "HCLUST Test with palindrom, any .bam file and .gff file Optimum clusters"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 -o results/HCLUST/SIMPLE_OPTIMUM
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 --nopal -o results/HCLUST/SIMPLE_OPTIMUM_nopal
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 --bam samples/mapping/reads.bam -o results/HCLUST/SIMPLE_OPTIMUM_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 --nopal --bam samples/mapping/reads.bam -o results/HCLUST/SIMPLE_OPTIMUM_nopal_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/SIMPLE_OPTIMUM_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/SIMPLE_OPTIMUM_nopal_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/SIMPLE_OPTIMUM_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa  -c 10 --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/SIMPLE_OPTIMUM_nopal_bam_cd

./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 -o results/HCLUST/F1000_OPTIMUM
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 --nopal -o results/HCLUST/F1000_OPTIMUM_nopal
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 --bam samples/mapping/reads.bam -o results/HCLUST/F1000_OPTIMUM_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 --nopal --bam samples/mapping/reads.bam -o results/HCLUST/F1000_OPTIMUM_nopal_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_OPTIMUM_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_OPTIMUM_nopal_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_OPTIMUM_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -c 10 --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_OPTIMUM_nopal_bam_cd

./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 -o results/HCLUST/PCA95_OPTIMUM
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 --nopal -o results/HCLUST/PCA95_OPTIMUM_nopal
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 --bam samples/mapping/reads.bam -o results/HCLUST/PCA95_OPTIMUM_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 --nopal --bam samples/mapping/reads.bam -o results/HCLUST/PCA95_OPTIMUM_nopal_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/PCA95_OPTIMUM_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/PCA95_OPTIMUM_nopal_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/PCA95_OPTIMUM_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa --pca 0.95 -c 10 --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/PCA95_OPTIMUM_nopal_bam_cd

echo "HCLUST Test with palindrom, any .bam file and .gff file"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall -o results/HCLUST/F1000_ball_hall
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn -o results/HCLUST/F1000_dunn
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette -o results/HCLUST/F1000_silhouette
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index -o results/HCLUST/F1000_ch-index
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index -o results/HCLUST/F1000_db-index

echo "HCLUST Test without palindrom, any .bam file and .gff file"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal -o results/HCLUST/F1000_ball_hall_nopal
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal -o results/HCLUST/F1000_dunn_nopal
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal -o results/HCLUST/F1000_silhouette_nopal
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal -o results/HCLUST/F1000_ch-index_nopal
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal -o results/HCLUST/F1000_db-index_nopal

echo "HCLUST Test with palindrom, any .gff file, with .bam file"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --bam samples/mapping/reads.bam -o results/HCLUST/F1000_ball_hall_bam
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --bam samples/mapping/reads.bam -o results/HCLUST/F1000_dunn_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --bam samples/mapping/reads.bam -o results/HCLUST/F1000_silhouette_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --bam samples/mapping/reads.bam -o results/HCLUST/F1000_ch-index_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --bam samples/mapping/reads.bam -o results/HCLUST/F1000_db-index_bam

echo "HCLUST Test without palindrom, any .gff file, with .bam file"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --bam samples/mapping/reads.bam -o results/HCLUST/F1000_ball_hall_nopal_bam
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --bam samples/mapping/reads.bam -o results/HCLUST/F1000_dunn_nopal_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --bam samples/mapping/reads.bam -o results/HCLUST/F1000_silhouette_nopal_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --bam samples/mapping/reads.bam  -o results/HCLUST/F1000_ch-index_nopal_bam
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --bam samples/mapping/reads.bam -o results/HCLUST/F1000_db-index_nopal_bam

echo "HCLUST Test with palindrom, any .bam file, with .gff file"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ball_hall_cd
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_dunn_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_silhouette_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ch-index_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_db-index_cd

echo "HCLUST Test without palindrom, any .bam file, with .gff file"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ball_hall_nopal_cd
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_dunn_nopal_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_silhouette_nopal_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ch-index_nopal_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_db-index_nopal_cd


echo "HCLUST Test with palindrom, with .bam and .gff files"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ball_hall_bam_cd
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_dunn_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_silhouette_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ch-index_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_db-index_bam_cd

echo "HCLUST Test without palindrom, with .bam and .gff files"
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ball-hall --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ball_hall_nopal_bam_cd
# ./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m dunn --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_dunn_nopal_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m silhouette --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_silhouette_nopal_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m ch-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_ch-index_nopal_bam_cd
./scripts/exploration/hierarchical_clustering.py -i samples/contigs/final.contigs.fa -f 1000 -m db-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/HCLUST/F1000_db-index_nopal_bam_cd
