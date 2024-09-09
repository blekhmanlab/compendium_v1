# for fast.adonis
library(plyr)
library(parallel)
# the rest
library(dplyr)
library(vegan)
library(fast.adonis)
library(rhdf5)

ITERS_PER_ADONIS <- 500
CPUS <- 1

start <- date()
print('Loading')

samplenames <- read.delim('keeper_names_sorted.txt')
full <- h5read('filtered.weighted.dm.pc','matrix')
rownames(full) <- samplenames$srs

metadata.raw <- read.delim('../tech.txt')

metadata <- metadata.raw[match(rownames(full), metadata.raw$srs),]

print('adonis')
asdf <- fast.adonis(full ~ worldregion+amplicon+avg_len+beating,
                    data=metadata, parallel=CPUS, permutations=ITERS_PER_ADONIS)

print(asdf[['aov.tab']])

results <- data.frame(
  factor = c('worldregion','amplicon','avg_len','beating'),
  SumOfSqs = asdf[['aov.tab']]$SumsOfSqs[1:4],
  R2 = asdf[['aov.tab']]$R2[1:4],
  F = asdf[['aov.tab']][["F.Model"]][1:4],
  P = asdf[['aov.tab']][["Pr(>F)"]][1:4]
)
done <- date()
print(start)
print(done)
write.csv(results, file='results.csv',
          row.names=F, quote=F)
