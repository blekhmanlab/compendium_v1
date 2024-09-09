# NOTE: This relies on a variable, points.annotated, that is
# defined and manipulated in visualization/figure2.R. It's
# easier to run this script after that one to avoid bouncing
# back and forth.

library(clusterSim)
library(dplyr)
library(magrittr) # for mod()

working <- points.annotated %>%
  filter(!region=='unknown') %>%
  filter(!region=='Oceania') %>%
  dplyr::select(mds1, mds2, mds3, mds4, mds5, mds6, mds7, mds8, region) %>%
  mutate(regionnum=as.integer(region))

table(working$regionnum)

regionlabels <- as.integer(working$regionnum)

working <- working %>%
  dplyr::select(mds1, mds2, mds3, mds4, mds5, mds6, mds7, mds8)

real <- index.DB(working, regionlabels, centrotypes="centroids")

results = c()
set.seed(711)
date()
for(i in seq(0,250000)) {
  if(mod(i,10000) == 0) {
    print(date())
    print(i)
  }

  newlabels <- sample(regionlabels)
  runtest <- index.DB(working, newlabels, centrotypes="centroids")
  results <- append(results, runtest$DB)
}
date()

results <- as.data.frame(results)

doneplot <- ggplot(results, aes(x=results)) +
  geom_histogram(bins=75) +
  geom_vline(xintercept=real$DB, color='red', linewidth=1) +
  labs(x='Davies-Bouldin Index', y='Iterations') +
  theme_bw()
ggsave('db_bootstrap.pdf', plot=doneplot, device='pdf',
       height=7, width=9, units='in')

# calculate percentile
findp <- ecdf(results$results)
findp(real$DB)
min(results$results)
