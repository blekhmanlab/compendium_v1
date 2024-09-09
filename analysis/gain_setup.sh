#!/bin/bash

# Make region-specific tables for subsampling
qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Eastern and South-Eastern Asia"' \
    --o-filtered-table regional_data/ese_asia.even1k.qza

qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Sub-Saharan Africa"' \
    --o-filtered-table regional_data/africa.even1k.qza

qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Europe and Northern America"' \
    --o-filtered-table regional_data/europe.even1k.qza

qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Northern Africa and Western Asia"' \
    --o-filtered-table regional_data/w_asia.even1k.qza

qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Australia/New Zealand"' \
    --o-filtered-table regional_data/australia.even1k.qza

qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Latin America and the Caribbean"' \
    --o-filtered-table regional_data/latam.even1k.qza

qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Central and Southern Asia"' \
    --o-filtered-table regional_data/c_asia.even1k.qza

qiime feature-table filter-samples \
    --i-table raw/even1k.biom.qza \
    --m-metadata-file raw/sample_metadata.tsv \
    --p-where '[worldregion]="Oceania"' \
    --o-filtered-table regional_data/oceania.even1k.qza

qiime feature-table group --i-table regional_data/no_unknowns.even1k.qza --m-metadata-file raw/sample_metadata.tsv --p-axis 'sample' --p-mode 'sum' --m-metadata-column 'worldregion' --o-grouped-table regional_data/grouped.even1k.qza

# make sure it worked
#qiime feature-table summarize --i-table grouped.even1k.qza --o-visualization grouped.even1k.qzv

# https://forum.qiime2.org/t/introducing-greengenes2-2022-10/25291
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.phylogeny.id.nwk.qza
mv 2022.10/2022.10.phylogeny.id.nwk.qza raw/

# https://github.com/biocore/unifrac
qiime diversity alpha-phylogenetic --i-table regional_data/grouped.even1k.qza \
--i-phylogeny raw/2022.10.phylogeny.id.nwk.qza \
--o-alpha-diversity base_diversity.qza \
--p-metric faith_pd --verbose

qiime tools export --input-path base_diversity.qza --output-path base_diversity
mv base_diversity/alpha-diversity.tsv results/base.tsv
