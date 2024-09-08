getPhylum <- function(taxon.Name) {
  temp <- unlist(strsplit(taxon.Name,split="[.]"))[2]
  if (identical(temp,integer(0))) {
    return(taxon.Name)
  }
  return(temp) #otherwise return the phylum
}
#############################
# SETTING UP THE DATA
#############################


topTaxa <- c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota")

phylum_scale <- setNames(c('#BCD2EE','#832161','#06D6A0','#E88873','#6153CC','gray'),
                         c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota","other") )

x <- readRDS(file="filtered.rds")

#read in metadata
metadata <- readRDS('metadata_from_rpackage.rds')

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
    plot.margin = unit(c(15,0,5,15),"pt")
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
stacked_bar <- panel_a + regions + plot_layout(ncol=1,nrow=2,heights=c(39,1))
