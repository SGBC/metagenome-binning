#! /usr/bin/env bash

echo "AF Test with palindrom, any .bam file and .gff file"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall -o results/AF/F2000_ball_hall
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn -o results/AF/F2000_dunn
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette -o results/AF/F2000_silhouette
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index -o results/AF/F2000_ch-index
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index -o results/AF/F2000_db-index

echo "AF Test without palindrom, any .bam file and .gff file"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal -o results/AF/F2000_ball_hall_nopal
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal -o results/AF/F2000_dunn_nopal
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal -o results/AF/F2000_silhouette_nopal
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal -o results/AF/F2000_ch-index_nopal
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal -o results/AF/F2000_db-index_nopal

echo "AF Test with palindrom, any .gff file, with .bam file"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --bam samples/mapping/reads.bam -o results/AF/F2000_ball_hall_bam
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --bam samples/mapping/reads.bam -o results/AF/F2000_dunn_bam
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --bam samples/mapping/reads.bam -o results/AF/F2000_silhouette_bam
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --bam samples/mapping/reads.bam -o results/AF/F2000_ch-index_bam
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --bam samples/mapping/reads.bam -o results/AF/F2000_db-index_bam

echo "AF Test without palindrom, any .gff file, with .bam file"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal --bam samples/mapping/reads.bam -o results/AF/F2000_ball_hall_nopal_bam
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal --bam samples/mapping/reads.bam -o results/AF/F2000_dunn_nopal_bam
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal --bam samples/mapping/reads.bam -o results/AF/F2000_silhouette_nopal_bam
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal --bam samples/mapping/reads.bam  -o results/AF/F2000_ch-index_nopal_bam
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal --bam samples/mapping/reads.bam -o results/AF/F2000_db-index_nopal_bam

echo "AF Test with palindrom, any .bam file, with .gff file"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ball_hall_cd
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_dunn_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_silhouette_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ch-index_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_db-index_cd

echo "AF Test without palindrom, any .bam file, with .gff file"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ball_hall_nopal_cd
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_dunn_nopal_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_silhouette_nopal_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ch-index_nopal_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_db-index_nopal_cd


echo "AF Test with palindrom, with .bam and .gff files"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ball_hall_bam_cd
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_dunn_bam_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_silhouette_bam_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ch-index_bam_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_db-index_bam_cd

echo "AF Test without palindrom, with .bam and .gff files"
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ball-hall --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ball_hall_nopal_bam_cd
# ./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m dunn --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_dunn_nopal_bam_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m silhouette --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_silhouette_nopal_bam_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m ch-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_ch-index_nopal_bam_cd
./scripts/exploration/affinity_propagation.py -i samples/contigs/final.contigs.fa -f 2000 -m db-index --nopal --bam samples/mapping/reads.bam --cd samples/prot_map/contigs_genes.gff -o results/AF/F2000_db-index_nopal_bam_cd
