# Human Microbiome Compendium v1

## Code supplement

This repository contains code related to the manuscript "Integration of 168,000 samples reveals global patterns of the human gut microbiome" by Abdill, Graham et al.

### Processing
* The `/pipeline` directory contains code used for retrieving, processing and consolidating the raw data from the Sequence Read Archive.
* `/visualization/setup.R` contains helper functions used in the scripts below.
* `/analysis/make_filtered_data.R` requires one input file, `taxonomic_table.csv`. This can be downloaded [from Zenodo](https://doi.org/10.5281/zenodo.8186993) (as `taxonomic_table.csv.gz`), decompressed, and used without modification. This script generates the data files used in all the other scripts.

### Analysis
The `/analysis` directory contains code used to generate and evaluate data for the project.
* `pcoa.R` contains the code used for the principal coordinates analysis in Figure 2.
* `rarefaction.R` contains the code used for the taxonomic discovery rate analysis in Figure 1I.
* `cluster_evaluation.R` contains the code used for the bootstrap analysis of clustering strength described in the manuscript.
* `pca.R` contains the code used for the principal components analysis used to determine regional signatures described in the manuscript. It relies on one external file, `sample_metadata.tsv`, available in the paper's [associated Zenodo repository](https://doi.org/10.5281/zenodo.8186993).
* `country_inference_check.R` contains the code used for the manual evaluation of the accuracy of the world region inference steps. The power calculation is first, followed by the procedure used to generate the randomly selected samples to validate.

### Visualization
The `/visualization` directory contains the R code used to generate the figures in our manuscript.

* `map_setup.R` lists the steps for installing the dependencies for generating the map in Figure 2A.
* `figure1.R` generates the panels in Figure 1 and associated supplementary material. It requires one external file, `rarefaction.rds`, that is stored in the `/data` directory.
* `visualization/figure2.R` generates the panels in Figure 2 and associated supplementary material. It requires several external files:
  * In the `data/` directory:
      * `rarefaction_diversity.rds`
      * `regions.csv`
  * From the paper's [associated Zenodo repository](https://doi.org/10.5281/zenodo.8186993):
    * `sample_metadata.tsv`
  * Generated bo `pcoa.R`:
    * `nmds.rds`
    * `pcoa_points.rds`
* `compendium_functions.R` sets up and imports libraries used for figures 3 and 4.
* `figure3.R` generates the panels in Figure 3 and associated supplementary material. It requires one external file, `unfiltered_rarefaction_by_read.rds`, that is stored in the `/data` directory.
* `figure4.R` generates the panels in Figure 4. It requires several external files, all available in the `/data` directory:
  * `diff_taxa_counts_for_4A.rds`
  * `fig4A_labels.rds`
  * `diff_abundant_pvalues_for_4B.rds`
  * `metadata_for_diffAbundance20230919.rds`
  * `taxon_names.tsv`

--
If you have any questions, please contact corresponding author Ran Blekhman at blekhman (at) uchicago.edu. Thanks.
