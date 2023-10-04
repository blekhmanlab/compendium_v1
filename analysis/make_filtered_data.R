library(dplyr)
library(data.table)
library(tidyr) # pivot_longer

MIN_READS_SAMPLE <- 10000
MIN_READS_TAXON <- 1000
MIN_SAMPLES_TAXON <- 100
MAX_WEIRD_PROPORTION <- 0.1
MAX_ARCHAEA_PROPORTION <- 0.1

# Read the data file
raw <- read.csv('taxonomic_table.csv')

# Set the row names to be sample names
x <- raw
rownames(x) <- x$sample
# Format the row names to be ${SAMPLE}_${PROJECT}
rownames(x) <- gsub('(\\w+)_consolidated.tsv(_\\w+)$', '\\1\\2', rownames(x))
############
x <- select(x, !c('X','sample'))

nrow(x) # total samples to start
ncol(x) # total taxa to start
nrow(x[rowSums(x) < MIN_READS_SAMPLE,]) # how many samples don't have enough reads?
x <- x[rowSums(x) >= MIN_READS_SAMPLE,] # get rid of them


# taxa must have >= 1000 reads across ALL samples
# (this is redundant to the next step, but is useful
# for thinning out the data frame before we start summing
# up columns)
ncol(x[,colSums(x) < MIN_READS_TAXON]) # taxa w less than 1000 reads total
x <- x[,colSums(x) >= MIN_READS_TAXON]

# taxa must appear in at least 100 samples
prev <- x %>% mutate_if(is.numeric, ~1 * (. > 0))
ncol(x[,colSums(prev) < MIN_SAMPLES_TAXON])
x <- x[,colSums(prev) >= MIN_SAMPLES_TAXON]

# after we get rid of taxa, make sure all our sample still have some reads left
nrow(x[rowSums(x) < MIN_READS_SAMPLE,])
x <- x[rowSums(x) >= MIN_READS_SAMPLE,]
rm(prev)

# Figure out the phylum of each taxon, and remove
# any for which the phylum is "NA"
taxphylum <- x
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
taxphylum <- combine_taxa(taxphylum) %>% make_rel()
taxphylum$weird <- taxphylum$NA.NA + taxphylum$Eukaryota.NA +
  taxphylum$Bacteria.NA + taxphylum$Archaea.NA
nrow(taxphylum[taxphylum$weird > MAX_WEIRD_PROPORTION,]) # samples to remove
x <- x[taxphylum$weird <= MAX_WEIRD_PROPORTION,]
rm(taxphylum)
x.unique <- unique(x)


######
# NOTE This is the ONLY data file with the
#      duplicate samples removed
#######
saveRDS(x.unique, file='filtered.rds')

taxphylum <- x
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
taxphylum <- combine_taxa(taxphylum)

write.csv(taxphylum, "phylum_complete.csv", row.names=T)
saveRDS(taxphylum, file='phylum.rds')
rm(taxphylum)
# -----------------
taxclass <- x
colnames(taxclass) <- gsub('^(\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxclass))
taxclass <- combine_taxa(taxclass)

write.csv(taxclass, "class_complete.csv", row.names=T)
saveRDS(taxclass, file='class.rds')
rm(taxclass)
# -----------------
taxorder <- x
colnames(taxorder) <- gsub('^(\\w+\\.\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxorder))
taxorder <- combine_taxa(taxorder)

saveRDS(taxorder, file='order.rds')
rm(taxorder)
# -----------------
taxfamily <- x
colnames(taxfamily) <- gsub('^(\\w+\\.\\w+\\.\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxfamily))
taxfamily <- combine_taxa(taxfamily)

saveRDS(taxfamily, file='family.rds')
rm(taxphylum)
# -----------------
taxgenus <- x
taxgenus <- combine_taxa(taxgenus)

saveRDS(taxgenus, file='2023/genus.rds')
rm(taxgenus)
