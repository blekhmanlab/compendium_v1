############ THIS SCRIPT IS FOR COMPARING DIVERSITY BETWEEN REGIONS

library(dplyr)
library(ggplot2)
library(patchwork)
library(vegan)

working <- readRDS('filtered.rds')

meta <- read.delim('sample_metadata.tsv', sep='\t')
meta$join <- paste(meta$project, meta$srr, sep='_')

totals <- as.data.frame(rowSums(working))
totals$sample <- rownames(totals)
totals <- rename(totals, total='rowSums(working)')

data <- meta %>%
  select(join, region) %>%
  inner_join(totals, by=join_by(join==sample)) %>%
  filter(region != 'unknown') %>%
  filter(region != '')

# make a version of the dataset with only samples that
# have a regional annotation

working$sample <- rownames(working)
only_regional <- working %>%
  inner_join(meta, by=join_by(sample==join)) %>%
  filter(region != 'unknown') %>%
  filter(region != '') %>%
  select(!any_of(colnames(meta)))
rownames(only_regional) <- only_regional$sample
only_regional$sample <- NULL


only_regional.annotated <- only_regional %>%
  mutate(join=rownames(only_regional)) %>%
  inner_join(data, by='join')
rownames(only_regional.annotated) <- only_regional.annotated$join

unrarefied <- only_regional.annotated %>%
  select(!c('join','region','total')) %>%
  diversity(index='shannon') %>%
  as.data.frame() %>%
  rename(shannon='.') %>%
  mutate(region=only_regional.annotated$region)

rawdiv <- ggplot(unrarefied, aes(x=reorder(region, shannon), y=shannon)) +
  geom_violin(draw_quantiles = 0.5) +
  coord_flip() +
  labs(y='Shannon diversity', x='Region') +
  theme_bw()

results <- data.frame(
  region=character(),
  diversity=double()
)
date()
set.seed(1234)
for(i in seq(1000)) {
  if(i %% 50 == 0) {
    print(i)
  }
  round <- only_regional.annotated %>%
    group_by(region) %>%
    slice_sample(n=1000) %>%
    ungroup()
  samplenames <- round$join
  
  round <- round %>%
    select(!c('join','region','total')) %>%
    rrarefy(10000) %>%
    diversity(index='shannon') %>%
    as.data.frame() %>%
    rename(shannon='.')
  round$sample <- samplenames
  
  summarized <- inner_join(round, meta, by=join_by(sample==join)) %>%
    group_by(region) %>%
    summarise(diversity=mean(shannon))
  
  results <- rbind(results, summarized)
}
date()

saveRDS(results, 'rarefaction_diversity.rds')

ggplot(results, aes(x=region, y=diversity)) +
  #geom_violin(draw_quantiles=c(0.5)) +
  geom_jitter(height=0) +
  coord_flip() +
  theme_bw()

ggplot(results, aes(x=reorder(region, diversity), y=diversity)) +
  #geom_violin(draw_quantiles=c(0.5)) +
  geom_violin(draw_quantiles = 0.5) +
  coord_flip() +
  labs(y='Shannon diversity', x='Region') +
  theme_bw()

overall <- results %>%
  group_by(region) %>%
  summarise(mean=mean(diversity), sd=sd(diversity))
