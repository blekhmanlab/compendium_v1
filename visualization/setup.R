library(dplyr)

# Given a table of relative abundances, create
# a new table in which each row shows a single
# taxon and how many samples it appears in
get_prevalence <- function(taxtable) {
  richness <- taxtable %>% mutate_if(is.numeric, ~1 * (. > 0))
  prevalence <- data.frame(colnames(richness), colSums(richness, na.rm=T))
  colnames(prevalence) <- c('taxon','samples')
  rownames(prevalence) <- NULL
  return(prevalence)
}

# From a count table in which each row is a sample,
# the input to this function is a single row of the
# table. It returns a row of the same size, but each
# count is converted into a decimal.
reldivide = function(x){
  x / sum(x)
}

# Convert an entire count table into a relative
# abundance table
make_rel <- function(table) {
  return(data.frame(t(apply(table, 1, reldivide))))
}

# Given a list of taxon names from DADA2, reformat
# them to look nice
format_names <- function(names) {
  # turn periods into spaces
  names <- gsub('.+\\.([^\\.]+)$', '\\1', names)
  names <- gsub('NA', '(Unassigned)', names)
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
