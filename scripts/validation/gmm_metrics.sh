#! /usr/bin/env bash

# ./scripts/validation/fast_metrics.py -i results/GMM/Optimal_f1500/* -m samples/seq_index.json --map -i graphs/VBGMM/ -o graphs/GMM/ -n "VBGMM_Optimal_f1500"

echo "GMM Test with palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic" 
echo "GMM Test without palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall_nopal/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall_nopal" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn_nopal/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn_nopal" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette_nopal/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette_nopal" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index_nopal" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index_nopal" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic_nopal/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic_nopal" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic_nopal/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic_nopal" 

echo "GMM Test with palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall_bam" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic_bam"
echo "GMM Test without palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall_nopal_bam" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic_nopal_bam"

echo "GMM Test with palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall_cd" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic_cd"
echo "GMM Test without palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall_nopal_cd" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic_nopal_cd" 

echo "GMM Test with palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette _bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic_bam_cd"
echo "GMM Test without palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ball_hall_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ball_hall_nopal_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/GMM/F2000_dunn_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_dunn_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_silhouette_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_silhouette_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_ch-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_ch-index_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_db-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_db-index_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_aic_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_aic_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/GMM/F2000_bic_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/GMM/ -n  "GMM_F2000_bic_nopal_bam_cd" 
