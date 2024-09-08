library(dplyr)
library(data.table)
library(tidyr) # pivot_longer

MIN_READS_SAMPLE <- 10000
MIN_READS_TAXON <- 1000
MIN_SAMPLES_TAXON <- 100
MAX_WEIRD_PROPORTION <- 0.1
MAX_ARCHAEA_PROPORTION <- 0.1

# When multiple columns have the same name,
# combine those columns by adding together
# the numbers in each
combine_taxa <- function(dataset) {
  dataset$srr <- rownames(dataset)
  dt <- dataset %>%
    pivot_longer(!srr, names_to = "taxon", values_to = "count") %>%
    data.table()
  final <- dt[, lapply(.SD, sum), by=list(srr, taxon)] %>%
    pivot_wider(names_from=taxon, values_from=count, values_fill=0) %>%
    data.frame() # tibbles don't like row names
  rownames(final) <- final$srr
  final$srr <- NULL
  return(final)
}

# Given a table of relative abundances, create
# a new table in which each row shows a single
# taxon and how many samples it appears in
get_prevalence <- function(taxtable) {
  richness <- taxtable %>% mutate_if(is.numeric, ~1 * (. > 0))
  prevalence <- data.frame(colnames(richness), colSums(richness))
  colnames(prevalence) <- c('taxon','samples')
  rownames(prevalence) <- NULL
  #prevalence$prop <- prevalence$samples / nrow(taxtable)
  return(prevalence)
}

# From a count table in which each row is a sample,
# the input to this function is a single row of the
# table. It returns a row of the same size, but each
# count is converted into a decimal.
rel = function(x){
  x / sum(x)
}

# Convert an entire count table into a relative
# abundance table
make_rel <- function(table) {
  return(data.frame(t(apply(table, 1, rel))))
}

# Similar to get_prevalence, this function takes
# a count table and returns a table in which each
# row lists a single taxon with its average abundance.
get_abundance <- function(taxtable) {
  taxtable.rel <- as.data.frame(t(apply(taxtable, 1, rel)))
  abundance <- data.frame(colnames(taxtable.rel), colSums(taxtable.rel, na.rm=TRUE)/nrow(taxtable.rel))
  rownames(abundance) <- NULL
  colnames(abundance) <- c('taxon','mean_abundance')
  return(abundance)
}

write.csv(raw, file='taxonomic_table.csv') # This is the taxonomic table we share with people


# Set the row names to be sample names
x <- as.data.frame(raw)
# row names are ${PROJECT}_${SAMPLE}
rownames(x) <- x$sample
x$sample <- NULL
x$X <- NULL
saveRDS(x, file='unfiltered.rds')

############
nrow(x) # total samples to start
ncol(x) # total taxa to start
nrow(x[rowSums(x) < MIN_READS_SAMPLE,]) # how many samples don't have enough reads?
x <- x[rowSums(x) >= MIN_READS_SAMPLE,] # get rid of empty samples


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
# Double-check that removing those samples didn't
# send any taxa below the thresholds
ncol(x[,colSums(x) < MIN_READS_TAXON]) # taxa w less than 1000 reads total
prev <- x %>% mutate_if(is.numeric, ~1 * (. > 0))
ncol(x[,colSums(prev) < MIN_SAMPLES_TAXON]) # taxa in less than 100 samples
rm(prev)

# Figure out the phylum of each taxon, and remove
# any for which the phylum is "NA"
taxphylum <- x
#colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
colnames(taxphylum) <- gsub('^([^\\.]+\\.[^\\.]+).+$', '\\1', colnames(taxphylum))

taxphylum <- combine_taxa(taxphylum) %>% make_rel()
taxphylum$weird <- taxphylum$NA.NA + taxphylum$Eukaryota.NA +
  taxphylum$Bacteria.NA + taxphylum$Archaea.NA
nrow(taxphylum[taxphylum$weird > MAX_WEIRD_PROPORTION,]) # samples to remove
x <- x[taxphylum$weird <= MAX_WEIRD_PROPORTION,]


######
# NOTE This is the ONLY data file with the
#      duplicate samples removed
#######
saveRDS(x, file='filtered.rds')

taxphylum <- x
colnames(taxphylum) <- gsub('^([^\\.]+\\.[^\\.]+).+$', '\\1', colnames(taxphylum))

taxphylum <- combine_taxa(taxphylum)
saveRDS(taxphylum, file='phylum.rds')
rm(taxphylum)
# -----------------
taxclass <- x
colnames(taxclass) <- gsub('^(\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxclass))
taxclass <- combine_taxa(taxclass)
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
rm(taxfamily)
