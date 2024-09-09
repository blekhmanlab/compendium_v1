import evident
import pandas as pd

metadata = pd.read_table("alpha-diversity.tsv", sep="\t", index_col=0)

with open('keeper_names_sorted.txt','r') as infile:
    infile.readline() # skip header
    keepers = [x[:-1] for x in infile]

faith_pd = metadata.query("index in @keepers")["faith_pd"]

metadata = pd.read_table("tech.txt", sep="\t", index_col=0)
adh = evident.UnivariateDataHandler(faith_pd, metadata, max_levels_per_category=15)

adh.calculate_effect_size(column="amplicon")
adh.calculate_effect_size(column="worldregion")
adh.calculate_effect_size(column="beating")
