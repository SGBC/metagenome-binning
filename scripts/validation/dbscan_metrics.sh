#! /usr/bin/env bash

echo "DBSCAN Test with palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index"

echo "DBSCAN Test without palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall_nopal/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall_nopal" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn_nopal/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn_nopal" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette_nopal/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette_nopal" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index_nopal" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index_nopal"

echo "DBSCAN Test with palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall_bam" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn_bam" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette_bam" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index_bam" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index_bam"

echo "DBSCAN Test without palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall_nopal_bam" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index_nopal_bam"


echo "DBSCAN Test with palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall_cd" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index_cd"

echo "DBSCAN Test without palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall_nopal_cd" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index_nopal_cd"


echo "DBSCAN Test with palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn_bam_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette_bam_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index_bam_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index_bam_cd"

echo "DBSCAN Test without palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ball_hall_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ball_hall_nopal_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_dunn_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_dunn_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_silhouette_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_silhouette_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_ch-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_ch-index_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/DBSCAN/F2000_db-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/DBSCAN/ -n  "DBSCAN_F2000_db-index_nopal_bam_cd"

