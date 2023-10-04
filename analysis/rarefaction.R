library(dplyr)
library(data.table)
library(tidyr) # pivot_longer

sessionInfo()

# ----------- RAREFACTION

REPLICATES <- 100
todo <- c(1, 10, 100, 500, 1000,
          seq(from=2000, to=10000, by=2000),
          seq(from=15000,to=35000,by=5000),
          seq(from=50000,to=110000,by=15000),
          125000, 150000
)

rarefaction_sampling <- function(taxtable) {
  richness <- taxtable
  richness <- richness %>% mutate_if(is.numeric, ~1 * (. > 0))
  found <- data.frame(colSums(richness)) %>% mutate_if(is.numeric, ~1 * (. > 0))
  found$taxon <- rownames(found)
  colnames(found) <- c('observed','taxon')
  observed <- sum(found$observed)
  rarefaction<-data.frame(nrow(richness), observed)
  names(rarefaction)<-c("scount","observed")

  for(scount in todo) {
    print(paste('SAMPLE SIZE:', scount))
    for(iter in seq(0,REPLICATES)) {
      sub <- dplyr::sample_n(richness, scount)
      found <- data.frame(colSums(sub)) %>%
        mutate_if(is.numeric, ~1 * (. > 0)) %>%
        rename(observed=colSums.sub.) %>%
        mutate(taxon=rownames(found))

      observed <- sum(found$observed)
      result <- data.frame(scount, observed)
      names(result)<-c("scount","observed")
      rarefaction <- rbind(rarefaction, result)
    }
  }

  print(paste('SAMPLE SIZE: ALL'))
  sub <- richness
  found <- data.frame(colSums(sub)) %>%
    mutate_if(is.numeric, ~1 * (. > 0)) %>%
    rename(observed=colSums.sub.) %>%
    mutate(taxon=rownames(found))

  observed <- sum(found$observed)
  result <- data.frame(nrow(richness), observed)
  names(result)<-c("scount","observed")
  rarefaction <- rbind(rarefaction, result)

  return(rarefaction)
}

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

raredata <- readRDS('unfiltered.rds')

# CONSOLIDATE COLUMN NAMES AT DIFFERENT TAXONOMIC LEVELS:
taxphylum <- raredata
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\..+$', '\\1', colnames(taxphylum))
taxphylum <- combine_taxa(taxphylum)
length(unique(colnames(taxphylum)))

taxclass <- raredata
colnames(taxclass) <- gsub('^(\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxclass))
taxclass <- combine_taxa(taxclass)

taxorder <- raredata
colnames(taxorder) <- gsub('^(\\w+\\.\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxorder))
taxorder <- combine_taxa(taxorder)

taxfamily <- raredata
colnames(taxfamily) <- gsub('^(\\w+\\.\\w+\\.\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxfamily))
taxfamily <- combine_taxa(taxfamily)
#---------------------------------------


set.seed(42)
print('phylum')
rarefaction.phylum <- rarefaction_sampling(taxphylum)
rarefaction.phylum$level <- 'phylum'
print('class')
rarefaction.class <- rarefaction_sampling(taxclass)
rarefaction.class$level <- 'class'
print('order')
rarefaction.order <- rarefaction_sampling(taxorder)
rarefaction.order$level <- 'order'
print('family')
rarefaction.family <- rarefaction_sampling(taxfamily)
rarefaction.family$level <- 'family'
print('genus')
rarefaction.genus <- rarefaction_sampling(raredata)
rarefaction.genus$level <- 'genus'

print('compiling')
rarefaction <- rbind(rarefaction.phylum, rarefaction.class,
                     rarefaction.order, rarefaction.family,
                     rarefaction.genus)
rarefaction$level <- as.factor(rarefaction$level)

rm(rarefaction.phylum)
rm(rarefaction.class)
rm(rarefaction.order)
rm(rarefaction.family)
rm(rarefaction.genus)

print('saving')
saveRDS(rarefaction, file='rarefaction.rds')
