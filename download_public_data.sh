#!/usr/bin/env bash
set -euo pipefail
mkdir -p data
cd data
echo "Downloading GSE342984 filtered 10x H5 files"
curl -L -C - -o GSM9944363_Sample1_filtered_feature_bc_matrix.h5 https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9944nnn/GSM9944363/suppl/GSM9944363_Sample1_filtered_feature_bc_matrix.h5
curl -L -C - -o GSM9944364_Sample2_filtered_feature_bc_matrix.h5 https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9944nnn/GSM9944364/suppl/GSM9944364_Sample2_filtered_feature_bc_matrix.h5
curl -L -C - -o GSM9944365_Sample3_filtered_feature_bc_matrix.h5 https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9944nnn/GSM9944365/suppl/GSM9944365_Sample3_filtered_feature_bc_matrix.h5
curl -L -C - -o GSM9944366_Sample4_filtered_feature_bc_matrix.h5 https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9944nnn/GSM9944366/suppl/GSM9944366_Sample4_filtered_feature_bc_matrix.h5
echo "Sample1/2 = 11 mM control; Sample3/4 = 33 mM high glucose"
