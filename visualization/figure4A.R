library(data.table)
library(tidyverse)
library(compositions)
library(RColorBrewer)
library(patchwork)
library(TreeSummarizedExperiment)

theme_set(theme_bw(base_size = 11))

#make a dataframe that contains the region names and their corresponding abbreviations
abbrs <- c("ANZ","CSA","ESEA","ENA","LAC","NAWA","SSA")
regions <- c("Australia/New Zealand","Central and Southern Asia",
             "Eastern and South-Eastern Asia","Europe and Northern America",
             "Latin America and the Caribbean", "Northern Africa and Western Asia",
             "Sub-Saharan Africa")
regionDF <<- data.frame(cbind(abbrs,regions))

region_scale <- setNames(c("#000000", "#e69d00", "#56b3e9", "#009e74",
                           "#f0e442","#d55e00", "#0071b2", "#cc79a7"),
                         c(
                           'Oceania',
                           'Northern Africa and Western Asia',
                           'Central and Southern Asia',
                           'Australia/New Zealand',
                           'Latin America and the Caribbean',
                           'Sub-Saharan Africa',
                           'Eastern and South-Eastern Asia',
                           'Europe and Northern America'
                         ))

phylum_scale <- setNames(c('#BCD2EE','#832161','#06D6A0','#E88873','#6153CC','gray'),
                         c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota","other") )

# Set up panels B and C
source('figure4B.R')
source('figure4C.R')


#read in the results of the discovery rate analysis
#this file is a data frame with four columns: sample_size, unique_taxa, total_reads, region
#unique taxa is the number of distinct taxa we found in n random microbiome samples from the given region for a single iteration
#total reads is the number of reads in all of the selected microbiome samples
res <- readRDS('unfiltered_rarefaction_by_read.rds')


res$region_sample <- paste(res$region,res$sample_size,sep="_")
remove <- which(is.na(res$unique_taxa)) #get rid of any sample size/region combos we don't have results for
res <- res[-remove,]


res$millions <- res$total_reads/1000000 #calculate the number of reads we have in millions of reads
res$taxaPerRead <- res$unique_taxa/res$total_reads #calculate the number of taxa we found per each read
res$taxaPerMillion <- res$unique_taxa/res$millions #calculate the number of taxa we found per one million reads


for (i in unique(res$region_sample)) {
  currRows <- which(res$region_sample == i)
  currData <- res[currRows,]
  currRegion <- unique(currData$region)
  currSampleSize <- unique(currData$sample_size)
  mean.reads <- mean(currData$total_reads)
  mean.taxa <- mean(currData$unique_taxa)
  taxaPerMillion <- mean(currData$taxaPerMillion)
  taxaPerRead <- mean(currData$taxaPerRead)
  if (i == unique(res$region_sample)[1]) {
    output <- data.frame(cbind(currRegion,currSampleSize,mean.reads,mean.taxa,taxaPerMillion,taxaPerRead))
  } else {
    temp <- data.frame(cbind(currRegion,currSampleSize,mean.reads,mean.taxa,taxaPerMillion,taxaPerRead))
    output <- rbind(output,temp)
  }
}

#each of these stats need to be numeric for plotting
output$currSampleSize <- as.numeric(output$currSampleSize)
output$mean.taxa <- as.numeric(output$mean.taxa)
output$mean.reads <- as.numeric(output$mean.reads)
output$taxaPerMillion <- as.numeric(output$taxaPerMillion)
output$taxaPerRead <- as.numeric(output$taxaPerRead)

#create a data frame to store the information for where each region label will go on the plot
textPos <- output
#these positions were decided through trial and error and are specific to our plotting area 
for (i in unique(textPos$currRegion)) {
  currMax <- max(textPos$currSampleSize[textPos$currRegion == i])
  remove <- which((textPos$currRegion == i) & (textPos$currSampleSize < currMax))
  textPos <- textPos[-remove,]
  if (i == "Australia/New Zealand") {
    textPos$mean.taxa[textPos$currRegion == i] <- textPos$mean.taxa[textPos$currRegion == i] + 100
  } else if (i == "Northern Africa and Western Asia") {
    textPos$mean.taxa[textPos$currRegion == i] <- textPos$mean.taxa[textPos$currRegion == i] - 50
  }
}

#create the inset plot in figure 4A
inset <- ggplot() + 
  geom_point(data=output,aes(y=taxaPerMillion,x=currSampleSize, color=currRegion)) + 
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) +
  scale_x_continuous(labels=scales::comma_format()) + 
  scale_y_continuous(labels=scales::comma_format()) + 
  geom_line(data=output,aes(y=taxaPerMillion,x=currSampleSize, color=currRegion)) +
  scale_y_log10()+
  xlab("Sample Size") +
  ylab("Taxa per Million Reads")  +
  theme(legend.position="none",
        axis.title = element_text(size=8)) 

#create the main plot in figure 4A
full <- ggplot() + 
  geom_point(data=output,aes(y=mean.taxa,x=currSampleSize,color=currRegion)) +
  geom_line(data=output,aes(y=mean.taxa,x=currSampleSize,color=currRegion)) +
  scale_x_continuous(labels=scales::comma_format()) + 
  scale_y_continuous(labels=scales::comma_format()) + 
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) +
  theme(legend.position = "none") + ylab("Unique taxa") + xlab("Sample size")# +

#combine the inset and the main panel
discovery_curve <- full + inset_element(inset,0.4,0.05,.98,0.65,align_to = "panel",ignore_tag = T)

l3 <- "
AAAABB
AAAABB
AAAABB
CCCCCC
CCCCCC
DDDDDD
DDDDDD
EEEEEE
EEEEEE
EEEEEE
FFFFFF
"
p3 <- discovery_curve + colError1 + phylum_plots + stacked_bar + networks + network.leg + plot_layout(design=l3)
