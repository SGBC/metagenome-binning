#! /usr/bin/env bash

# echo "HCLUST Test with palindrom, any .bam file and .gff file Optimum clusters"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM/* -n "HCLUST_SIMPLE_OPTIMUM"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM_nopal/* -n "HCLUST_SIMPLE_OPTIMUM_nopal"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM_bam/* -n "HCLUST_SIMPLE_OPTIMUM_bam"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM_nopal_bam/* -n "HCLUST_SIMPLE_OPTIMUM_nopal_bam"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM_cd/* -n "HCLUST_SIMPLE_OPTIMUM_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM_nopal_cd/* -n "HCLUST_SIMPLE_OPTIMUM_nopal_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM_bam_cd/* -n "HCLUST_SIMPLE_OPTIMUM_bam_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/SIMPLE_OPTIMUM_nopal_bam_cd/* -n "HCLUST_SIMPLE_OPTIMUM_nopal_bam_cd"

# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM/* -n "HCLUST_F2000_OPTIMUM"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM_nopal/* -n "HCLUST_F2000_OPTIMUM_nopal"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM_bam/* -n "HCLUST_F2000_OPTIMUM_bam"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM_nopal_bam/* -n "HCLUST_F2000_OPTIMUM_nopal_bam"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM_cd/* -n "HCLUST_F2000_OPTIMUM_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM_nopal_cd/* -n "HCLUST_F2000_OPTIMUM_nopal_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM_bam_cd/* -n "HCLUST_F2000_OPTIMUM_bam_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_OPTIMUM_nopal_bam_cd/* -n "HCLUST_F2000_OPTIMUM_nopal_bam_cd"

# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM/* -n "HCLUST_PCA95_OPTIMUM"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM_nopal/* -n "HCLUST_PCA95_OPTIMUM_nopal"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM_bam/* -n "HCLUST_PCA95_OPTIMUM_bam"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM_nopal_bam/* -n "HCLUST_PCA95_OPTIMUM_nopal_bam"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM_cd/* -n "HCLUST_PCA95_OPTIMUM_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM_nopal_cd/* -n "HCLUST_PCA95_OPTIMUM_nopal_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM_bam_cd/* -n "HCLUST_PCA95_OPTIMUM_bam_cd"
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/PCA95_OPTIMUM_nopal_bam_cd/* -n "HCLUST_PCA95_OPTIMUM_nopal_bam_cd"

echo "HCLUST Test with palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall/* -n "HCLUST_F2000_ball_hall" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn/* -n "HCLUST_F2000_dunn" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette/* -n "HCLUST_F2000_silhouette" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index/* -n "HCLUST_F2000_ch-index" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index/* -n "HCLUST_F2000_db-index"

echo "HCLUST Test without palindrom, any .bam file and .gff file"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall_nopal/* -n "HCLUST_F2000_ball_hall_nopal" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn_nopal/* -n "HCLUST_F2000_dunn_nopal" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette_nopal/* -n "HCLUST_F2000_silhouette_nopal" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index_nopal/* -n "HCLUST_F2000_ch-index_nopal" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index_nopal/* -n "HCLUST_F2000_db-index_nopal"

echo "HCLUST Test with palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall_bam/* -n "HCLUST_F2000_ball_hall_bam" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn_bam/* -n "HCLUST_F2000_dunn_bam" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette_bam/* -n "HCLUST_F2000_silhouette_bam" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index_bam/* -n "HCLUST_F2000_ch-index_bam" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index_bam/* -n "HCLUST_F2000_db-index_bam"

echo "HCLUST Test without palindrom, any .gff file, with .bam file"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall_nopal_bam/* -n "HCLUST_F2000_ball_hall_nopal_bam" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn_nopal_bam/* -n "HCLUST_F2000_dunn_nopal_bam" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette_nopal_bam/* -n "HCLUST_F2000_silhouette_nopal_bam" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index_nopal_bam/* -n "HCLUST_F2000_ch-index_nopal_bam" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index_nopal_bam/* -n "HCLUST_F2000_db-index_nopal_bam"

echo "HCLUST Test with palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall_cd/* -n "HCLUST_F2000_ball_hall_cd" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn_cd/* -n "HCLUST_F2000_dunn_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette_cd/* -n "HCLUST_F2000_silhouette_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index_cd/* -n "HCLUST_F2000_ch-index_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index_cd/* -n "HCLUST_F2000_db-index_cd"

echo "HCLUST Test without palindrom, any .bam file, with .gff file"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall_nopal_cd/* -n "HCLUST_F2000_ball_hall_nopal_cd" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn_nopal_cd/* -n "HCLUST_F2000_dunn_nopal_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette_nopal_cd/* -n "HCLUST_F2000_silhouette_nopal_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index_nopal_cd/* -n "HCLUST_F2000_ch-index_nopal_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index_nopal_cd/* -n "HCLUST_F2000_db-index_nopal_cd"


echo "HCLUST Test with palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall_bam_cd/* -n "HCLUST_F2000_ball_hall_bam_cd" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn_bam_cd/* -n "HCLUST_F2000_dunn_bam_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette_bam_cd/* -n "HCLUST_F2000_silhouette_bam_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index_bam_cd/* -n "HCLUST_F2000_ch-index_bam_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index_bam_cd/* -n "HCLUST_F2000_db-index_bam_cd"

echo "HCLUST Test without palindrom, with .bam and .gff files"
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ball_hall_nopal_bam_cd/* -n "HCLUST_F2000_ball_hall_nopal_bam_cd" &
# ./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_dunn_nopal_bam_cd/* -n "HCLUST_F2000_dunn_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_silhouette_nopal_bam_cd/* -n "HCLUST_F2000_silhouette_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_ch-index_nopal_bam_cd/* -n "HCLUST_F2000_ch-index_nopal_bam_cd" &
./scripts/validation/fast_metrics.py -m samples/seq_index.json --map --nopal -o graphs/HCLUST/ -i results/HCLUST/F2000_db-index_nopal_bam_cd/* -n "HCLUST_F2000_db-index_nopal_bam_cd"
