#!/bin/sh
set -ex

iter=$1

mkdir iterations/run${iter}
cd iterations/run${iter}

# BUILD SUBSAMPLED ASV TABLES
mkdir subsamples
qiime feature-table subsample-ids --i-table ../../regional_data/ese_asia.even1k.qza \
    --p-subsampling-depth 1000 \
    --p-axis sample \
    --o-sampled-table subsamples/subsample.ese_asia.even1k.qza
qiime feature-table subsample-ids --i-table ../../regional_data/africa.even1k.qza \
    --p-subsampling-depth 1000 \
    --p-axis sample \
    --o-sampled-table subsamples/subsample.africa.even1k.qza
qiime feature-table subsample-ids --i-table ../../regional_data/europe.even1k.qza \
    --p-subsampling-depth 1000 \
    --p-axis sample \
    --o-sampled-table subsamples/subsample.europe.even1k.qza
qiime feature-table subsample-ids --i-table ../../regional_data/w_asia.even1k.qza \
    --p-subsampling-depth 1000 \
    --p-axis sample \
    --o-sampled-table subsamples/subsample.w_asia.even1k.qza
qiime feature-table subsample-ids --i-table ../../regional_data/australia.even1k.qza \
    --p-subsampling-depth 1000 \
    --p-axis sample \
    --o-sampled-table subsamples/subsample.australia.even1k.qza
qiime feature-table subsample-ids --i-table ../../regional_data/latam.even1k.qza \
    --p-subsampling-depth 1000 \
    --p-axis sample \
    --o-sampled-table subsamples/subsample.latam.even1k.qza
qiime feature-table subsample-ids --i-table ../../regional_data/c_asia.even1k.qza \
    --p-subsampling-depth 1000 \
    --p-axis sample \
    --o-sampled-table subsamples/subsample.c_asia.even1k.qza

qiime feature-table merge \
--i-tables subsamples/subsample.ese_asia.even1k.qza \
    subsamples/subsample.africa.even1k.qza \
    subsamples/subsample.europe.even1k.qza \
    subsamples/subsample.w_asia.even1k.qza \
    subsamples/subsample.australia.even1k.qza \
    subsamples/subsample.latam.even1k.qza \
    subsamples/subsample.c_asia.even1k.qza \
--p-overlap-method error_on_overlapping_sample \
--o-merged-table subsample.even1k.qza

rm -rf subsamples

qiime feature-table group --i-table subsample.even1k.qza \
    --m-metadata-file ../../raw/sample_metadata.tsv --p-axis 'sample' --p-mode 'sum' \
    --m-metadata-column 'worldregion' --o-grouped-table subsample.grouped.even1k.qza

# Build the tables that combine REGIONS together
qiime feature-table group --i-table subsample.grouped.even1k.qza --m-metadata-file ../../combinations.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'rep1' --o-grouped-table rep1.even1k.qza
qiime feature-table group --i-table subsample.grouped.even1k.qza --m-metadata-file ../../combinations.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'rep2' --o-grouped-table rep2.even1k.qza
qiime feature-table group --i-table subsample.grouped.even1k.qza --m-metadata-file ../../combinations.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'rep3' --o-grouped-table rep3.even1k.qza
qiime feature-table group --i-table subsample.grouped.even1k.qza --m-metadata-file ../../combinations.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'rep4' --o-grouped-table rep4.even1k.qza
qiime feature-table group --i-table subsample.grouped.even1k.qza --m-metadata-file ../../combinations.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'rep5' --o-grouped-table rep5.even1k.qza
qiime feature-table group --i-table subsample.grouped.even1k.qza --m-metadata-file ../../combinations.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'rep6' --o-grouped-table rep6.even1k.qza
qiime feature-table group --i-table subsample.grouped.even1k.qza --m-metadata-file ../../combinations.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'rep7' --o-grouped-table rep7.even1k.qza

qiime feature-table merge \
--i-tables subsample.grouped.even1k.qza \
    rep1.even1k.qza \
    rep2.even1k.qza \
    rep3.even1k.qza \
    rep4.even1k.qza \
    rep5.even1k.qza \
    rep6.even1k.qza \
    rep7.even1k.qza \
--p-overlap-method error_on_overlapping_sample \
--o-merged-table subsample.region_combos.even1k.qza

# CALCULATE BASE DIVERSITY
qiime diversity alpha-phylogenetic --i-table subsample.region_combos.even1k.qza \
--i-phylogeny ../../raw/2022.10.phylogeny.id.nwk.qza \
--o-alpha-diversity diversity.qza \
--p-metric faith_pd

# GET RESULTS OUT
qiime tools export --input-path diversity.qza --output-path interim
mv interim/alpha-diversity.tsv ../../results/${iter}.tsv
rm -rf interim
