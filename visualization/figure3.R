################
# Figure 3B
################

# this script requires filtered.rds and sample_metadata.tsv
source('compendium_functions.R')
x <- readRDS('filtered.rds')

# Figure out the phylum of each taxon, and remove
# any for which the phylum is "NA"
taxphylum <- x
#a few taxon names have dashes or spaces that mess up the parsing, so remove those
colnames(taxphylum) <- gsub("-","",colnames(taxphylum))
colnames(taxphylum) <- gsub(" ","",colnames(taxphylum))

#shorten the taxon name to just the kingdom.phylum
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
data <- combine_taxa(taxphylum) 

# Find the five most prevalent phyla so we can focus on those
top_phyla <- get_prevalence(data)
top_phyla$taxon <- gsub('^\\w+\\.(\\w+)$', '\\1', top_phyla$taxon)
colnames(data) <- gsub('^\\w+\\.(\\w+)$', '\\1', colnames(data))

top_phyla <- top_phyla %>%
  arrange(desc(samples)) %>%
  slice_head(n=5) %>%
  rbind(data.frame(
    taxon=c('other'),
    samples=c(0)
  ))

# Start with the phylum-level data, and turn
# it into a relative abundance table
final.rel <- make_rel(data)

# find the taxa with the highest variance, so we can use them
toplot <- final.rel %>%
  mutate(sample=rownames(final.rel)) %>%
  pivot_longer(!sample, names_to = "taxon", values_to = "rel")

# If it's not a top phylum, don't include it as a name
toplot[!(toplot$taxon %in% top_phyla$taxon[1:5]),]$taxon <- 'other'

###### REGION-LEVEL HISTOGRAMS

meta <- read.delim('sample_metadata.tsv', sep='\t') %>%
  mutate(sample=paste(project, srr, sep='_'))

toplot.sum <- toplot %>% left_join(meta,by="sample")
keep <- toplot.sum$region %in% regionDF$regions
toplot.sum <- toplot.sum[keep,]

toplot.sum$region <- gsub("South-Eastern","S.E.",toplot.sum$region)
toplot.sum$region <- gsub("Eastern","E.",toplot.sum$region)
toplot.sum$region <- gsub("Southern", "S.",toplot.sum$region)
toplot.sum$region <- gsub("Northern","N.",toplot.sum$region)
toplot.sum$region <- gsub("Western","W.",toplot.sum$region)
toplot.sum$region <- gsub("Latin America","Latin Amer.",toplot.sum$region)


toplot.sum$taxon <- gsub("Bacteroidota","Bacteroidetes",toplot.sum$taxon)
toplot.sum$taxon <- gsub("Actinobacteriota","Actinobacteria",toplot.sum$taxon)
names(phylum_scale) <-gsub("Actinobacteriota","Actinobacteria",names(phylum_scale))
names(phylum_scale) <-gsub("Bacteroidota","Bacteroidetes",names(phylum_scale))

phylum_plots <- ggplot(toplot.sum, aes(x=rel, y=after_stat(count), color=taxon)) +
  geom_freqpoly(size=1) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),limits=c(1,5000)) + #if I set the scale min at zero it can't plot
  scale_x_continuous(labels = c("",scales::percent(.5),"",scales::percent(1)), limits=c(0, 1.0), expand=c(0,0),breaks=c(.25,.5,.75,1)) +
  theme(plot.margin = unit(c(20,0,0,0),"pt"),
        axis.title = element_text(size=11),
        legend.position="bottom") +
  labs(
    x='Relative abundance',
    y='Samples'
  ) + 
  facet_wrap(~region,nrow=1, labeller = labeller(region = label_wrap_gen(13))) + 
  scale_fill_manual(values=phylum_scale, aesthetics=c('colour','fill')) + 
  guides(color=guide_legend(title="Phylum"))

################
# Figure 3C
################
getPhylum <- function(taxon.Name) {
  temp <- unlist(strsplit(taxon.Name,split="[.]"))[2]
  if (identical(temp,integer(0))) {
    return(taxon.Name)
  }
  return(temp) #otherwise return the phylum
}

topTaxa <- c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota")

phylum_scale <- setNames(c('#BCD2EE','#832161','#06D6A0','#E88873','#6153CC','gray'),
                         c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota","other") )

x <- readRDS(file="filtered.rds")

#read in metadata
metadata <- read.delim('sample_metadata.tsv', sep='\t') %>%
  mutate(sample=paste(project, srr, sep='_'))

#get rid of unnecessary samples and reorder metadata so it matches the samples in x
metadata <- metadata[metadata$sample %in% rownames(x),]
ord <- match(rownames(x),metadata$sample)
metadata <- metadata[ord,]

#get rid of any samples without a region assignment
keep <- metadata$region %in% regionDF$regions
metadata <- metadata[keep,]
x <- x[keep,]

#collapse x to the phylum level
taxphylum <- x
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
colnames(taxphylum) <- sapply(colnames(taxphylum),getPhylum)
final.rel <- combine_taxa(taxphylum) %>% make_rel()

final.rel$region <- metadata$region

if (exists('myData')) {
  rm(myData)
}
for (i in 1:nrow(regionDF)) {
  currRegionSamples <- which(final.rel$region == regionDF$regions[i]) 
  subset <- final.rel[currRegionSamples,] #get samples from current region of interest
  set.seed(42)
  subset <- subset[sample(nrow(subset),round(0.1*nrow(subset))),] #take 10% of those samples
  
  #we only needed the region info to split up into the correct samples
  #but we can get rid of it now
  subset <- subset(subset,select=-c(region))
  
  # if we want to order the samples using certain taxa,
  # we need an extra copy of those, which we'll attach
  # to every entry for every sample in the pivot_longer version.
  subset$sample <- rownames(subset)
  
  
  # now go BACK to long form to plot
  final.long <- subset %>% pivot_longer(!sample, names_to = "taxon", values_to = "rel")
  
  final.long[!(final.long$taxon %in% topTaxa),]$taxon <- 'other'
  final.long$taxon <- factor(final.long$taxon, levels=c(topTaxa, 'other'))
  
  #order the samples by top taxa
  final.long$sample <- factor(final.long$sample, levels=subset[
    order(subset[[topTaxa[1]]],
          subset[[topTaxa[1]]]+subset[[topTaxa[2]]],
          subset[[topTaxa[1]]]+subset[[topTaxa[2]]]+subset[[topTaxa[3]]],
          subset[[topTaxa[1]]]+subset[[topTaxa[2]]]+subset[[topTaxa[3]]]+subset[[topTaxa[4]]],
          subset[[topTaxa[1]]]+subset[[topTaxa[2]]]+subset[[topTaxa[3]]]+subset[[topTaxa[4]]]+subset[[topTaxa[5]]]),
  ]$sample)
  
  final.long$region <- regionDF$regions[i]
  if (exists("myData")) {
    myData <- rbind(myData,final.long) #if the data frame exists, add this region's data to it
  } else {
    myData <- final.long #if the data frame doesn't exist, create it
  }
}


panel_a <- ggplot(myData, aes(fill=taxon, y=rel, x=sample)) +
  geom_bar(stat="identity",width=1) +
  theme_bw() +
  theme(
    legend.position='bottom',
    axis.text=element_text(size=11),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title = element_text(size=11),
    axis.title.x = element_blank(),
  ) +
  scale_fill_manual(values=phylum_scale, aesthetics=c('colour','fill'),
                    labels=c("Firmicutes","Proteobacteria","Actinobacteria","Bacteroidetes","Desulfobacterota","Other")) +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  labs(y='Relative abundance', fill='Phylum')

phylum_legend <- cowplot::get_legend(panel_a) 
panel_a <- panel_a + theme(legend.position = "none")

a_legend <- cowplot::get_legend(panel_a +
                                  theme(
                                    legend.position='bottom',
                                    legend.title = element_text(size=10),
                                    legend.text = element_text(size=7), #changed font size from 11 to 14 for presentation figure
                                    plot.tag.position  = c(1, 1)
                                  )
)


#we want to make a colored bar to go along the bottom x-axis with the region color labels
#to do that, we need to know the indices where each region starts in the main plot
#so iterate through the data frame we're plotting to find the first incidence of each region
breaks <- c(rep(NA,7))
for (i in 1:nrow(regionDF)) {
  currbreak <- min(which(myData$region == regionDF$region[i]))
  breaks[i] <- currbreak
}

#make the region label
regions <- ggplot(data=myData) + geom_tile(aes(x=1:nrow(myData),y=1,fill=region)) + scale_x_discrete(breaks=breaks,labels=regionDF$regions,expand=c(0,0)) + 
  theme(
    legend.position='none',
    axis.text=element_text(size=8),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    plot.title=element_text(hjust=0.95,vjust=-10),
    panel.spacing =  unit(c(0, 0, 0, 0), "pt"),
    plot.margin = unit(c(0,0,0,0),"pt")
  ) + scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) + xlab("World region") + scale_y_discrete(expand=c(0,0))

#add the region label as an inset to the main panel
stacked_bar <-  panel_a + inset_element(regions,0,-0.2,1,0,align_to = "panel",ignore_tag = TRUE) 


##############
# Figure 3A
#############
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

#create the inset plot in figure 3A
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

#create the main plot in figure 3A
full <- ggplot() + 
  geom_point(data=output,aes(y=mean.taxa,x=currSampleSize,color=currRegion)) +
  geom_line(data=output,aes(y=mean.taxa,x=currSampleSize,color=currRegion)) +
  scale_x_continuous(labels=scales::comma_format()) + 
  scale_y_continuous(labels=scales::comma_format()) + 
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) +
  theme(legend.position = "none") + ylab("Unique taxa") + xlab("Sample size") +
  geom_text(data=textPos[!(textPos$currRegion == "Europe and Northern America"),],
            aes(x=currSampleSize+1000,y=mean.taxa,color=currRegion,label=currRegion),hjust=0,size=3) +
  geom_text(data=textPos[(textPos$currRegion == "Europe and Northern America"),],
            aes(x=currSampleSize,y=mean.taxa+100,color=currRegion,label=currRegion),hjust=1,size=3) 

#combine the inset and the main panel
discovery_curve <- full + inset_element(inset,0.4,0.05,.98,0.65,align_to = "panel",ignore_tag = T)


#combine all three panels into a single plot
plot_panels <- discovery_curve + plot_spacer() + phylum_plots + stacked_bar + plot_layout(ncol=1,heights=c(5,-.3,2,2.5)) + 
  plot_annotation(tag_levels = "A",caption = "\n\n",theme=theme(plot.caption = element_text(hjust=0.5)))
