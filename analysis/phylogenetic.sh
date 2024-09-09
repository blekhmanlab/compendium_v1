#!/bin/bash -l
#SBATCH -J 168k-dataset
#SBATCH --time 24:00:00
#SBATCH --mem 300gb
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mail-type=FAIL,TIME_LIMIT_80,INVALID_DEPEND
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err

mamba activate qiime2-amplicon-2023.9
function logger () { 
    echo "$(date) :: ${*}"; 
    echo "$(date) :: ${*}" 1>&2; 
}

set -x 
set -e
set -o pipefail

# Faith PD results, unweighted and weighted UniFrac distances, their PCoA coordinates, and per-feature taxonomy data

#qiime greengenes2 non-v4-16s \
#    --i-table merged.w_md.biom.qza \
#    --i-sequences merged.w_md.fna.qza\
#    --i-backbone $HOME/greengenes2/release/2022.10/2022.10.backbone.full-length.fna.qza \
#    --o-mapped-table merged.w_md.gg2-2022.10-cref99.biom.qza \
#    --o-representatives merged.w_md.gg2-2022.10-cref99.reps.fna.qza \
#    --p-threads ${SLURM_CPUS_PER_TASK}
#
## alternatively, could use Naive Bayes mapping if desired
#qiime greengenes2 taxonomy-from-table \
#    --i-table merged.w_md.gg2-2022.10-cref99.biom.qza \
#    --i-reference-taxonomy $HOME/greengenes2/release/2022.10/2022.10.taxonomy.id.nwk.qza \
#    --o-classification merged.w_md.gg2-2022.10-cref99.tax.tsv.qza
#
#qiime feature-table rarefy \
#    --i-table merged.w_md.gg2-2022.10-cref99.biom.qza \
#    --p-sampling-depth 1000 \
#    --o-rarefied-table merged.w_md.gg2-2022.10-cref99.even1k.biom.qza

#qiime diversity alpha-phylogenetic \
#    --i-table merged.w_md.gg2-2022.10-cref99.even1k.biom.qza \
#    --i-phylogeny $HOME/greengenes2/release/2022.10/2022.10.phylogeny.id.nwk.qza \
#    --p-metric faith_pd \
#    --o-alpha-diversity merged.w_md.gg2-2022.10-cref99.even1k.faithpd.tsv.qza
#

# running via qiime2 is currently quite inefficient as it is not using the
# zero-copy semantics in the unifrac library. there are also deserialization
# costs to compute pcoa, which can be done natively in unifrac.
#qiime diversity beta-phylogenetic \
#    --i-table merged.w_md.gg2-2022.10-cref99.even1k.biom.qza \
#    --i-phylogeny $HOME/greengenes2/release/2022.10/2022.10.phylogeny.id.nwk.qza \
#    --p-metric unweighted_unifrac \
#    --o-distance-matrix merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.qza \
#    --p-threads ${SLURM_CPUS_PER_TASK}
#
#qiime diversity beta-phylogenetic \
#    --i-table merged.w_md.gg2-2022.10-cref99.even1k.biom.qza \
#    --i-phylogeny $HOME/greengenes2/release/2022.10/2022.10.phylogeny.id.nwk.qza \
#    --p-metric weighted_normalized_unifrac \
#    --o-distance-matrix merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.qza \
#    --p-threads ${SLURM_CPUS_PER_TASK}
#
#qiime diversity pcoa \
#    --i-distance-matrix merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.qza \
#    --o-pcoa merged.w_md.gg2-2022.10-cref99.even1k.unweighted.pc.qza \
#    --p-number-of-dimensions 10
#
#qiime diversity pcoa \
#    --i-distance-matrix merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.qza \
#    --o-pcoa merged.w_md.gg2-2022.10-cref99.even1k.weighted.pc.qza \
#    --p-number-of-dimensions 10
#qiime tools export \
#    --input-path merged.w_md.gg2-2022.10-cref99.even1k.biom.qza \
#    --output-path merged.w_md.gg2-2022.10-cref99.even1k 

#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#ssu -i merged.w_md.gg2-2022.10-cref99.even1k/feature-table.biom \
#    -t $HOME/greengenes2/release/2022.10/2022.10.phylogeny.id.nwk \
#    -m unweighted_fp32 \
#    -o $PANFS2/merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.pc \
#    --format hdf5_fp32 \
#    --pcoa 10 
#
#ssu -i merged.w_md.gg2-2022.10-cref99.even1k/feature-table.biom \
#    -t $HOME/greengenes2/release/2022.10/2022.10.phylogeny.id.nwk \
#    -m weighted_normalized_fp32 \
#    -o $PANFS2/merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.pc \
#    --format hdf5_fp32 \
#    --pcoa 10 
#
#cp $PANFS2/merged.w_md.gg2-2022.10-cref99.even1k.*.dm.pc .

#python export-pcoa-to-qza.py merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.pc &
#python export-pcoa-to-qza.py merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.pc &
#wait

#qiime emperor plot \
#    --i-pcoa merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.pc.qza \
#    --m-metadata-file merged.w_md.tsv \
#    --o-visualization merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.pc.qzv &
#
#qiime emperor plot \
#    --i-pcoa merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.pc.qza \
#    --m-metadata-file merged.w_md.tsv \
#    --o-visualization merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.pc.qzv &

qiime empress community-plot \
    --i-tree $HOME/greengenes2/release/2022.10/2022.10.phylogeny.id.nwk.qza \
    --i-feature-table merged.w_md.gg2-2022.10-cref99.even1k.biom.qza \
    --i-pcoa merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.pc.qza \
    --m-sample-metadata-file merged.w_md.tsv \
    --m-feature-metadata-file merged.w_md.gg2-2022.10-cref99.tax.tsv.qza \
    --o-visualization merged.w_md.gg2-2022.10-cref99.even1k.unweighted.dm.pc.empress.qzv &

qiime empress community-plot \
    --i-tree $HOME/greengenes2/release/2022.10/2022.10.phylogeny.id.nwk.qza \
    --i-feature-table merged.w_md.gg2-2022.10-cref99.even1k.biom.qza \
    --i-pcoa merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.pc.qza \
    --m-sample-metadata-file merged.w_md.tsv \
    --m-feature-metadata-file merged.w_md.gg2-2022.10-cref99.tax.tsv.qza \
    --o-visualization merged.w_md.gg2-2022.10-cref99.even1k.weighted.dm.pc.empress.qzv &
wait
