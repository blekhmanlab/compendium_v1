library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)

gm_mean = function(x){
  exp(mean(log(x[x > 0])))
}
rclr <- function(a) {
  answer <- log(a/gm_mean(a))
  #answer[is.infinite(answer)] <- 0
  answer[] <- lapply(answer, function(i) if(is.numeric(i)) ifelse(is.infinite(i), 0, i) else i)
  return(answer)
}
find_taxa <- function(regionselect, pc_cutoff, factor_cutoff) {
  print(date())
  set.seed(42)
  ohwow <- prcomp(
    data.transformed[regions==regionselect,],
    tol=0,
    center=FALSE
  )
  print(date())
  
  loadings <- as.data.frame(ohwow$rotation)
  
  var_explained = ohwow$sdev^2 /sum(ohwow$sdev^2)
  print('calculating scree data')
  scree <- data.frame(
    axis=ncol(data.transformed)+1-row_number(var_explained),
    explained=var_explained
  )
  scree$cumulative <- cumsum(scree$explained)
  
  # First, select only the PCs that explain more than 1 percent of the variation.
  # The evaluate the variance explained by each variable within those 10 PCs.
  print('Evaluating loadings')
  sums <- loadings %>%
    #dplyr::select(c(1:sum(var_explained > 0.01))) %>%
    dplyr::select(c(1:(1+sum(scree$cumulative < pc_cutoff)))) %>%
    mutate(explained=rowSums(.^2))  %>%
    arrange(explained)
  print(sum(sums$explained))
  return(sums)
}

working <- readRDS('filtered.rds')
working$join <- rownames(working)

meta <- read.delim('sample_metadata.tsv', sep='\t') %>%
  mutate(sample=paste(project, srr, sep='_'))

meta$join <- paste(meta$project, meta$srr, sep='_')

data <- meta %>%
  dplyr::select(join, region) %>%
  inner_join(working, by='join') %>%
  filter(!region == 'unknown') %>% 
  dplyr::select(!c('join')) %>%
  unique()

regions <- data$region
data <- data %>% dplyr::select(!c('region'))
rm(meta)
rm(working)

data.transformed <- rclr(data)
rm(data)

colnames(data.transformed) <- gsub('\\s','_', colnames(data.transformed))
colnames(data.transformed) <- gsub('\\[','_', colnames(data.transformed))
colnames(data.transformed) <- gsub('\\]','_', colnames(data.transformed))
colnames(data.transformed) <- gsub('-','_', colnames(data.transformed))
colnames(data.transformed) <- gsub('\\(','_', colnames(data.transformed))
colnames(data.transformed) <- gsub('\\)','_', colnames(data.transformed))




australia <- find_taxa('Australia/New Zealand', pc_cutoff=0.5, factor_cutoff=0.5)
csasia <- find_taxa('Central and Southern Asia', pc_cutoff=0.5, factor_cutoff=0.5)
ese_asia <- find_taxa('Eastern and South-Eastern Asia', pc_cutoff=0.5, factor_cutoff=0.5)
europe <- find_taxa('Europe and Northern America', pc_cutoff=0.5, factor_cutoff=0.5)
latam <- find_taxa('Latin America and the Caribbean', pc_cutoff=0.5, factor_cutoff=0.5)
nafrica <- find_taxa('Northern Africa and Western Asia', pc_cutoff=0.5, factor_cutoff=0.5)
safrica <- find_taxa('Sub-Saharan Africa', pc_cutoff=0.5, factor_cutoff=0.5)

get_topn <- function(data) {
  n = 10
  
  data %>%
    top_n(n, explained) %>%
    rownames() %>%
    return()
}

tokeep <- c(get_topn(australia), get_topn(csasia),
          get_topn(ese_asia), get_topn(europe),
          get_topn(latam), get_topn(nafrica),
          get_topn(safrica)
        ) %>%
        unique()
