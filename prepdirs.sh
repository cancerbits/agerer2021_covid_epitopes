#!/bin/bash
mkdir sars042
mkdir sars060
mkdir plots
cp ./GEOupload/barcodes042.tsv.gz ./sars042/barcodes042.tsv.gz
cp ./GEOupload/features042.tsv.gz ./sars042/features042.tsv.gz
cp ./GEOupload/matrix042.mtx.gz ./sars042/matrix042.mtx.gz
cp ./GEOupload/QC042rep.csv ./sars042/QC042rep.csv


cp ./GEOupload/barcodes060.tsv.gz ./sars060/barcodes060.tsv.gz
cp ./GEOupload/features060.tsv.gz ./sars060/features060.tsv.gz
cp ./GEOupload/matrix060.mtx.gz ./sars060/matrix060.mtx.gz
cp ./GEOupload/QC060rep.csv ./sars060/QC060rep.csv

mv ./sars060/barcodes060.tsv.gz ./sars060/barcodes.tsv.gz
mv ./sars060/features060.tsv.gz ./sars060/features.tsv.gz
mv ./sars060/matrix060.mtx.gz ./sars060/matrix.mtx.gz

mv ./sars042/barcodes042.tsv.gz ./sars042/barcodes.tsv.gz
mv ./sars042/features042.tsv.gz ./sars042/features.tsv.gz
mv ./sars042/matrix042.mtx.gz ./sars042/matrix.mtx.gz
