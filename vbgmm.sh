#! /usr/bin/env bash

echo "VBGMM + Filter 1500_OPTIMAL"

./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/Optimal_f1500 -c 10 > vbgmm.log
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/Optimal_f1500_bam -c 10 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/Optimal_f1500_bam_nopal -c 10 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/Optimal_f1500_cov -c 10 --cd samples/prot_map/contigs_prot.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --nopal -o results/VBGMM/Optimal_f1500_cov_nopal -c 10 --cd samples/prot_map/contigs_prot.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/Optimal_f1500_bam_cov -c 10 --cd samples/prot_map/contigs_prot.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/Optimal_f1500_bam_cov_nopal -c 10 --cd samples/prot_map/contigs_prot.gff 

echo "VBGMM + F1500 BIC"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/BIC_f1500 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/BIC_f1500_nopal --nopal 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_f1500_bam 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_f1500_bam_nopal 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/BIC_f1500_cd --cd samples/prot_map/contigs_prot.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --nopal -o results/VBGMM/BIC_f1500_cd_nopal --cd samples/prot_map/contigs_prot.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_f1500_bam_cd --cd samples/prot_map/contigs_prot.gff 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_f1500_bam_cd_nopal --cd samples/prot_map/contigs_prot.gff 

echo "VBGMM + F1500 AIC"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/AIC_f1500 --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/AIC_f1500_nopal --aic --nopal 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_f1500_bam --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_f1500_bam_nopal --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/AIC_f1500_cd --cd samples/prot_map/contigs_prot.gff --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --nopal -o results/VBGMM/AIC_f1500_cd_nopal --cd samples/prot_map/contigs_prot.gff --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_f1500_bam_cd --cd samples/prot_map/contigs_prot.gff --aic 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_f1500_bam_cd_nopal --cd samples/prot_map/contigs_prot.gff --aic 

echo "VBGMM + F1500 BIC"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/BIC_f1500_db --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/BIC_f1500_nopal_db --nopal --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_f1500_bam_db --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_f1500_bam_nopal_db --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/BIC_f1500_cd_db --cd samples/prot_map/contigs_prot.gff --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --nopal -o results/VBGMM/BIC_f1500_cd_nopal_db --cd samples/prot_map/contigs_prot.gff --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_f1500_bam_cd_db --cd samples/prot_map/contigs_prot.gff --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_f1500_bam_cd_nopal_db --cd samples/prot_map/contigs_prot.gff --dbscore 

echo "VBGMM + F1500 AIC DBSCORE"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/AIC_f1500_db --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/AIC_f1500_nopal_db --aic --nopal --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_f1500_bam_db --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_f1500_bam_nopal_db --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 -o results/VBGMM/AIC_f1500_cd_db --cd samples/prot_map/contigs_prot.gff --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --nopal -o results/VBGMM/AIC_f1500_cd_nopal_db --cd samples/prot_map/contigs_prot.gff --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/AIC_f1500_bam_cd_db --cd samples/prot_map/contigs_prot.gff --aic --dbscore 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/AIC_f1500_bam_cd_nopal_db --cd samples/prot_map/contigs_prot.gff --aic --dbscore 

echo "VBGMM + F1500 BIC BAM weight"
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_f1500_bam_weight -w 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_f1500_bam_nopal_weight -w 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam -o results/VBGMM/BIC_f1500_bam_cd_weight --cd samples/prot_map/contigs_prot.gff -w 
./scripts/exploration/vbgmm.py -i samples/contigs/final.contigs.fa -f 1500 --bam samples/mapping/reads.bam --nopal -o results/VBGMM/BIC_f1500_bam_cd_nopal_weight --cd samples/prot_map/contigs_prot.gff -w 

echo "done"