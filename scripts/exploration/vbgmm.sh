#! /usr/bin/env bash

echo "VBGMM + Filter 1500_OPTIMAL"

./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/Optimal_F1000 -c 10 > vbgmm.log
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/Optimal_F1000_bam -c 10 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/Optimal_F1000_bam_nopal -c 10 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/Optimal_F1000_cov -c 10 --cd samples/prot_map/contigs_genes.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --nopal -o results/VBGMM/Optimal_F1000_cov_nopal -c 10 --cd samples/prot_map/contigs_genes.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/Optimal_F1000_bam_cov -c 10 --cd samples/prot_map/contigs_genes.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/Optimal_F1000_bam_cov_nopal -c 10 --cd samples/prot_map/contigs_genes.gff 

echo "VBGMM + F1000 BIC"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/BIC_F1000 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/BIC_F1000_nopal --nopal 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_F1000_bam 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_F1000_bam_nopal 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/BIC_F1000_cd --cd samples/prot_map/contigs_genes.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --nopal -o results/VBGMM/BIC_F1000_cd_nopal --cd samples/prot_map/contigs_genes.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_F1000_bam_cd --cd samples/prot_map/contigs_genes.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_F1000_bam_cd_nopal --cd samples/prot_map/contigs_genes.gff 

echo "VBGMM + F1000 AIC"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/AIC_F1000 --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/AIC_F1000_nopal --aic --nopal 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_F1000_bam --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_F1000_bam_nopal --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/AIC_F1000_cd --cd samples/prot_map/contigs_genes.gff --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --nopal -o results/VBGMM/AIC_F1000_cd_nopal --cd samples/prot_map/contigs_genes.gff --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_F1000_bam_cd --cd samples/prot_map/contigs_genes.gff --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_F1000_bam_cd_nopal --cd samples/prot_map/contigs_genes.gff --aic 

echo "VBGMM + F1000 BIC"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/BIC_F1000_db --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/BIC_F1000_nopal_db --nopal --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_F1000_bam_db --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_F1000_bam_nopal_db --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/BIC_F1000_cd_db --cd samples/prot_map/contigs_genes.gff --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --nopal -o results/VBGMM/BIC_F1000_cd_nopal_db --cd samples/prot_map/contigs_genes.gff --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_F1000_bam_cd_db --cd samples/prot_map/contigs_genes.gff --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_F1000_bam_cd_nopal_db --cd samples/prot_map/contigs_genes.gff --dbscore 

echo "VBGMM + F1000 AIC DBSCORE"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/AIC_F1000_db --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/AIC_F1000_nopal_db --aic --nopal --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_F1000_bam_db --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_F1000_bam_nopal_db --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 -o results/VBGMM/AIC_F1000_cd_db --cd samples/prot_map/contigs_genes.gff --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --nopal -o results/VBGMM/AIC_F1000_cd_nopal_db --cd samples/prot_map/contigs_genes.gff --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_F1000_bam_cd_db --cd samples/prot_map/contigs_genes.gff --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_F1000_bam_cd_nopal_db --cd samples/prot_map/contigs_genes.gff --aic --dbscore 

echo "VBGMM + F1000 BIC BAM weight"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_F1000_bam_weight -w 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_F1000_bam_nopal_weight -w 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_F1000_bam_cd_weight --cd samples/prot_map/contigs_genes.gff -w 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1000 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_F1000_bam_cd_nopal_weight --cd samples/prot_map/contigs_genes.gff -w 

echo "done"