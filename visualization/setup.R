library(dplyr)

# Given a table of relative abundances, create
# a new table in which each row shows a single
# taxon and how many samples it appears in
get_prevalence <- function(taxtable) {
  richness <- taxtable %>% mutate_if(is.numeric, ~1 * (. > 0))
  prevalence <- data.frame(colnames(richness), colSums(richness, na.rm=T))
  colnames(prevalence) <- c('taxon','samples')
  rownames(prevalence) <- NULL
  #prevalence$prop <- prevalence$samples / nrow(taxtable)
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
