library(dplyr)

#### Sample size calculation

z <- 1.96 # 95 percent confidence interval
p <- 0.95 # what's our guess for the level of accuracy?
d <- 0.05 # precision level. Result would be (point estimate) plus or minus d
dropout <- 0.25 # what proportion will we need to throw out?

# SAMPLES REQUIRED
(((z^2) * p * (1-p)) / (d^2)) / (1-dropout)

############
##
# Big evaluation
##

meta <- read.delim('sample_metadata.tsv', sep='\t') %>%
  filter(region != 'unknown')
table(meta$region)
set.seed(123)
projects <- sample(unique(meta$project), 100)

set.seed(124)
to_eval <- meta %>%
  filter(project %in% projects) %>%
  group_by(project) %>%
  sample_n(size=1)

write.csv(to_eval, file='country_toeval.csv', row.names=F)
