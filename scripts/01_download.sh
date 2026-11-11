#!/usr/bin/env bash
# Download public, de-identified GEO processed matrices (no raw IDAT). Resumable (-C -).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw
B=https://ftp.ncbi.nlm.nih.gov/geo/series
get() { [ -s "data/raw/$2" ] && echo "have $2" || curl -fL -C - -o "data/raw/$2" "$1"; }
# GSE40279 (Hannum 2013): series matrix = sample metadata (age, gender, ...) + beta values (473k probes x 656 samples)
get $B/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz GSE40279_series_matrix.txt.gz
# GSE87571 (Johansson 2013): series matrix = metadata only; betas are in two supplementary matrices
get $B/GSE87nnn/GSE87571/matrix/GSE87571_series_matrix.txt.gz GSE87571_series_matrix.txt.gz
get $B/GSE87nnn/GSE87571/suppl/GSE87571_matrix1of2.txt.gz GSE87571_matrix1of2.txt.gz
get $B/GSE87nnn/GSE87571/suppl/GSE87571_matrix2of2.txt.gz GSE87571_matrix2of2.txt.gz
# depositor's additional sample characteristics (Sample_Group, year, smoking, genotype dosages, Houseman cell estimates, ...)
get $B/GSE87nnn/GSE87571/suppl/GSE87571_additional_sample_chararcteristics.xlsx GSE87571_additional_sample_chararcteristics.xlsx
# Reinius 2012 purified leukocytes (450K), used to build the cell-type reference (02b)
get $B/GSE35nnn/GSE35069/matrix/GSE35069_series_matrix.txt.gz GSE35069_series_matrix.txt.gz
(cd data/raw && gzip -t *.gz && ls -l && sha256sum *.gz *.xlsx > ../../results/raw_files.sha256)
