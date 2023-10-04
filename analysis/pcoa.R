library(dplyr)
library(data.table)
library(vegan) # for distance calc
library(bigmds) # https://arxiv.org/abs/2007.11919

sessionInfo()

# Given a table of relative abundances, create
# a new table in which each row shows a single
# taxon and how many samples it appears in
get_prevalence <- function(taxtable) {
  richness <- taxtable %>% mutate_if(is.numeric, ~1 * (. > 0))
  prevalence <- data.frame(colnames(richness), colSums(richness))
  colnames(prevalence) <- c('taxon','samples')
  rownames(prevalence) <- NULL
  return(prevalence)
}

gm_mean = function(x){
  exp(mean(log(x[x > 0])))
}
rclr <- function(a) {
    answer <- log(a/gm_mean(a))
    answer[is.infinite(answer)] <- 0
    return(answer)
}

print('Loading data')
reads <- readRDS('class.rds')

print('Removing duplicates')
reads <- unique(reads) # Required to prevent MDS from breaking

print('Determining prevalence')
prevalence.class <- get_prevalence(reads)

print('CALCULATING')
date()
set.seed(45)
# plotting 8 MDS dimensions:
nmds <- divide_conquer_mds(reads, 10000, 16, 8, n_cores = 16, dist_fn=vegdist, method='robust.aitchison')
date()
print('Saving results')
saveRDS(nmds, 'nmds.rds')
points <- data.frame(nmds$points)

print('Saving points')
points$sample <- rownames(reads)
colnames(points) <- c('mds1','mds2','mds3','mds4','mds5','mds6','mds7','mds8','sample')
points$project <- gsub('^(\\w+)_\\w+', '\\1', points$sample)

saveRDS(points, 'pcoa_points.rds')
