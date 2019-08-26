#! /usr/bin/env bash

# echo "AF Test with palindrom, any .bam file and .gff file Optimum clusters"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE_nopal/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE_nopal"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE_bam"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE_nopal_bam"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE_cd"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE_nopal_cd"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE_bam_cd"
# ./scripts/exploration/fast_metrics.py -i results/AF/SIMPLE_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_SIMPLE_nopal_bam_cd"


echo "AF Test with palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index"

echo "AF Test without palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall_nopal/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall_nopal" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn_nopal/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn_nopal"
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette_nopal/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette_nopal" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index_nopal" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index_nopal/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index_nopal"


echo "AF Test with palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall_bam" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn_bam" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette_bam" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index_bam" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index_bam"

echo "AF Test without palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall_nopal_bam" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index_nopal_bam" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index_nopal_bam/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index_nopal_bam"


echo "AF Test with palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall_cd" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index_cd"

echo "AF Test without palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall_nopal_cd" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn_nopal_cd"
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index_nopal_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index_nopal_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index_nopal_cd" 

echo "AF Test with palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn_bam_cd"
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette_bam_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index_bam_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index_bam_cd"

echo "AF Test without palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -i results/AF/F2000_ball_hall_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ball_hall_nopal_bam_cd" &
# ./scripts/validation/fast_metrics.py -i results/AF/F2000_dunn_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_dunn_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_silhouette_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_silhouette_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_ch-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_ch-index_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -i results/AF/F2000_db-index_nopal_bam_cd/* -m samples/seq_index.json --map --nopal -o graphs/AF/ -n "AF_F2000_db-index_nopal_bam_cd" 
