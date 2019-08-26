#! /usr/bin/env bash

# ./scripts/validation/fast_metrics.py -i results/VBGMM/Optimal_F2000/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_Optimal_F2000"
# ./scripts/validation/fast_metrics.py -i results/VBGMM/Optimal_F2000_bam/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_Optimal_F2000_bam"
# ./scripts/validation/fast_metrics.py -i results/VBGMM/Optimal_F2000_bam_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_Optimal_F2000_bam_nopal"
# ./scripts/validation/fast_metrics.py -i results/VBGMM/Optimal_F2000_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_Optimal_F2000_cd"
# ./scripts/validation/fast_metrics.py -i results/VBGMM/Optimal_F2000_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_Optimal_F2000_cd_nopal"
# ./scripts/validation/fast_metrics.py -i results/VBGMM/Optimal_F2000_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_Optimal_F2000_bam_cd"
# ./scripts/validation/fast_metrics.py -i results/VBGMM/Optimal_F2000_bam_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_Optimal_F2000_bam_cd_nopal"

./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_cd" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_cd_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_cd" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_cd_nopal" &

./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_cd" & 
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_cd_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam_cd" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam_cd_nopal" &

./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_db" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_nopal_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_nopal_db" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_db" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_nopal_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_nopal_db" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_cd_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_cd_db" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_cd_nopal_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_cd_nopal_db" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_cd_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_cd_db" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_cd_nopal_db/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_cd_nopal_db" &

./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_cd" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_cd_nopal" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam_cd" &
./scripts/validation/fast_metrics.py -i results/VBGMM/AIC_F2000_bam_cd_nopal/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_AIC_F2000_bam_cd_nopal" &

./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_weight/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_weight" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_nopal_weight/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_nopal_weight" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_cd_weight/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_cd_weight" &
./scripts/validation/fast_metrics.py -i results/VBGMM/BIC_F2000_bam_cd_nopal_weight/* -m samples/seq_index.json --map --nopal -o graphs/VBGMM/ -n "VBGMM_BIC_F2000_bam_cd_nopal_weight"
