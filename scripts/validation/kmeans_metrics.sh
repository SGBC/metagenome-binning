#! /usr/bin/env bash

# echo "K-MEANS Test with palindrom, any .bam file and .gff file Optimum clusters"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM_nopal/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM_nopal"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM_bam"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM_nopal_bam"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM_cd"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM_nopal_cd"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM_bam_cd"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/SIMPLE_OPTIMUM_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_SIMPLE_OPTIMUM_nopal_bam_cd"

# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM_nopal/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM_nopal"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM_bam"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM_nopal_bam"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM_cd"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM_nopal_cd"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM_bam_cd"
# ./scripts/exploration/fast_metrics.py -i results/K-MEANS/F2000_OPTIMUM_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_OPTIMUM_nopal_bam_cd"


echo "K-MEANS Test with palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index"

echo "K-MEANS Test without palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall_nopal/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall_nopal" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn_nopal/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn_nopal" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette_nopal/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette_nopal" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index_nopal" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index_nopal"


echo "K-MEANS Test with palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall_bam" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn_bam" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette_bam" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index_bam" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index_bam"

echo "K-MEANS Test without palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall_nopal_bam" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index_nopal_bam"


echo "K-MEANS Test with palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall_cd" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index_cd"

echo "K-MEANS Test without palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall_nopal_cd" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index_nopal_cd"

echo "K-MEANS Test with palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn_bam_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette_bam_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index_bam_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index_bam_cd"

echo "K-MEANS Test without palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ball_hall_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ball_hall_nopal_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_dunn_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_dunn_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_silhouette_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_silhouette_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_ch-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_ch-index_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/K-MEANS/F2000_db-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/K-MEANS/ -n "K-MEANS_F2000_db-index_nopal_bam_cd"
