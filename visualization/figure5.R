#load the input file
#this contains the number of taxa that are differentially abundant between each region-region pair
taxaCounts <- readRDS('diff_taxa_counts_for_5A.rds')
labels <- readRDS('fig4A_labels.rds')

plot_region_order <- c("Australia/New Zealand","Central and Southern Asia","Eastern and South-Eastern Asia","Latin America and the Caribbean",
                       "Northern Africa and Western Asia","Sub-Saharan Africa","Europe and Northern America")

NUM_TAXA <- 65 #total number of taxa included in the differential abundance analysis
bubble <- ggplot() + geom_point(data=taxaCounts,aes(x=factor(region1,levels=plot_region_order),y=factor(region2,level=plot_region_order),
                                                     size=freq,fill="grey"),shape=21,stroke=0) + 
  geom_text(data=taxaCounts,aes(x=factor(region1,levels=plot_region_order),y=factor(region2,levels=(plot_region_order)),label=freq),
            size=3,vjust=3) + 
  geom_label(data=labels,aes(x=factor(region2,levels=plot_region_order),y=factor(region2,levels=plot_region_order),label=regionLabel,fill=region2),size=2) +
  theme(legend.position = "top", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 9),
        panel.spacing = unit(c(0,0,0,0),"pt"),
        plot.margin = unit(c(0,0,0,0),"pt") ) +
  scale_radius(range = c(2, 13),
               limits=c(0,70),
               breaks=c(10,20,30,50,NUM_TAXA)) +
  guides(size = guide_legend(override.aes = list(fill = "grey", stroke = .25), 
                             label.position = "bottom",
                             title.position = "top", 
                             order = 1),
         fill=FALSE,color=FALSE) +
  labs(size = "Number of \nDifferentially \nAbundant Taxa",x = NULL, y = NULL) +
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill'))

#get the legend from the plot
size_legend <- cowplot::get_legend(bubble)

bubble <- bubble + theme(legend.position = 'none') + inset_element(size_legend,0.15,0.6,0.4,0.95,align_to="panel")

#create a label with the region colors to add to the x-axis of the plot
label <- ggplot(taxaCounts) + geom_tile(aes(x=factor(region2,levels=plot_region_order),y=1,fill=region2)) +  
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) + 
  theme(#axis.text.x=element_text(angle = 45, hjust=1,size=8),
    axis.text.x=element_blank(),
    legend.position="none",
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    panel.spacing = unit(c(0,0,0,0),"pt"),
    plot.margin = unit(c(0,0,0,0),"pt")
  ) + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))

#combine the blot with the region 
overview <- wrap_elements(bubble + label + plot_layout(heights = c(39,1)))


##INPUTS
#'diff_abundant_pvalues_for_5B.rds', which contains results of the differential abundance analysis
#'diff_abundant_taxa.rds' is a list of the taxa included in the differential abundance analysis
#'metadata_for_diffAbundance20230919.rds' is the metadata file used for the differential abundance analysis
#'taxon_names.tsv' contains a list of the full name of the taxon as listed in the compendium, and the parsed genus name

library(ggridges)
library(viridis)

pvals <- readRDS('diff_abundant_pvalues_for_5B.rds') #load the data containing the p values from the 
genus <- readRDS('genus.rds')
taxa_in_analysis <- readRDS('diff_abundant_taxa.rds') #list of taxa included in the differential abundance analysis
keepTaxa <- colnames(genus) %in% taxa_in_analysis #we want to get rid of taxa that were excluded from our analysis

metadata <- readRDS('metadata_for_diffAbundance.rds')
keepSamples <- rownames(genus) %in% metadata$sample #we want to get rid of samples that were excluded from our analysis
genus[genus == 0] <- 1
genus.rel <- make_rel(genus)

genus.rel <- genus.rel[keepSamples,keepTaxa]

#this file contains the parsed names of the taxa
taxaNames <- read.table(file="taxon_names.tsv",sep="\t",header=T)
taxaNames <- taxaNames[taxaNames$taxon %in% pvals$taxon,] #get rid of unnecessary taxa

#replace the full name with the parsed name in the pvals dataframe
for (i in 1:length(taxaNames$full.taxon)) {
  currFullTaxon <- taxaNames$full.taxon[i]
  currTaxon <- taxaNames$taxon[i]
  pvals$taxon[pvals$taxon == currFullTaxon] <- currTaxon
}


if (!(identical(metadata$sample,rownames(genus.rel)))) {
  stop("samples don't match")
}
ord <- match(colnames(genus.rel),taxaNames$full.taxon)
taxaNames <- taxaNames[ord,]
colnames(genus.rel) <- taxaNames$taxon #update the column names in genus to the parsed names

taxa.means <- apply(genus.rel,MARGIN=2,FUN=mean) #calculate the mean relative abundance of each taxon
taxa.means <- data.frame(taxa.means)
taxa.means$taxon <- rownames(taxa.means)

#sort the data frame with the average abundances in order from highest to lowest abundance
ord <- order(taxa.means$taxa.means,decreasing = T)
taxa.means <- taxa.means[ord,]

#create a barplot with the average relative abundances of each taxon
mean.bar <- ggplot(data=taxa.means) + geom_col(aes(x=taxa.means,y=factor(taxon,levels = rev(taxon)))) + xlab("Mean Relative \nAbundance") + 
  theme(axis.title.y=element_blank(),
        panel.spacing =  unit(c(0, 0, 0, 0), "pt"),
        plot.margin = unit(c(0,0,0,0),"pt"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x  = element_text(size=7),
        axis.text.x = element_text(size=6)) + scale_x_continuous(breaks=c(0,0.1))

#create a heat map with the p values comparing each region to Europe and Northern America for each taxon
pval.heatmap <- ggplot(data=pvals) + geom_tile(aes(x=region2,y=factor(taxon,levels=rev(taxa.means$taxon)),fill=log10(padj))) + 
  scale_fill_gradient2(low="red",high="white") + 
  theme(axis.text.x = element_blank(),
        axis.title=element_blank(), legend.position = "top",
        axis.ticks = element_blank(),
        text=element_text(size=8),legend.title = element_text(size=9)) + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0)) + 
  guides(fill=guide_colorbar(title="Adjusted P Value",title.position = "top",title.hjust = 0.5))

#create a label for the heatmap with the region color to identify it at the bottom
pval.label <- ggplot(data=pvals) + geom_tile(aes(x=region2,y=1,fill=region2)) + 
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) + theme(legend.position = "none",
                                                                                #  axis.text.x = element_text(angle = 45, vjust = 1,hjust=1),
                                                                                axis.text.x=element_blank(),
                                                                                axis.text.y = element_blank(),
                                                                                axis.ticks = element_blank(),
                                                                                axis.title = element_blank(),
                                                                                panel.spacing =  unit(c(0, 0, 0, 0), "pt"),
                                                                                plot.margin = unit(c(0,0,0,0),"pt")) +
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))


pval.heatmap.legend <- cowplot::get_legend(pval.heatmap)

pval <- pval.heatmap + pval.label + plot_layout(heights=c(24,1)) #combine the heatmap and the region label into one plot

genus.rel$region <- metadata$region
#now we want to caluclate the average relative abundance for each taxon in each region
for (i in 1:nrow(regionDF)) {
  currRegion <- regionDF$regions[i]
  counts <- genus.rel[genus.rel$region == currRegion,c(1:65)] #get a dataframe of the taxon abundances for only the samples in the current region
  means <- data.frame(apply(counts,MARGIN=2,mean)) %>% rownames_to_column(var="taxon") #calculate the mean for each taxon
  colnames(means) <- c("taxon","mean")
  means$region <- currRegion
  if (i == 1) {
    meanDF <- means #create the dataframe
  } else {
    meanDF <- rbind(meanDF,means) #if the dataframe already exists, add to it
  }
}

#make a heatmap showing the average relative abundance of each taxon in each region
mean.heatmap <- ggplot(meanDF) + geom_tile(aes(x=factor(region,levels = c("Australia/New Zealand","Central and Southern Asia",
                                                                          "Eastern and South-Eastern Asia","Latin America and the Caribbean",
                                                                          "Northern Africa and Western Asia","Sub-Saharan Africa",
                                                                          "Europe and Northern America")),
                                               y=factor(taxon,levels = rev(taxa.means$taxon)),
                                               fill=log10(mean))) + 
  scale_fill_viridis() +theme(
    axis.title.y=element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing =  unit(c(0, 0, 0, 0), "pt"),
    plot.margin = unit(c(0,0,0,0),"pt"),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9),
    legend.position="top") + geom_vline(xintercept=6.5,color="white",linetype="dashed") + 
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0)) + 
  guides(fill=guide_colorbar(title="Mean Relative Abundance",title.position = "top",title.hjust = 0.5))

mean.heatmap.legend <- cowplot::get_legend(mean.heatmap)

#make a region label for the x-axis
mean.region.label <- ggplot(meanDF) +  geom_tile(aes(x=factor(region,levels = c("Australia/New Zealand","Central and Southern Asia",
                                                                                "Eastern and South-Eastern Asia","Latin America and the Caribbean",
                                                                                "Northern Africa and Western Asia","Sub-Saharan Africa",
                                                                                "Europe and Northern America")),y=1, fill=region)) + 
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) + theme(legend.position = "none",
                                                                                axis.text.x=element_blank(),
                                                                                axis.text.y = element_blank(),
                                                                                axis.ticks = element_blank(),
                                                                                axis.title = element_blank(),
                                                                                panel.spacing =  unit(c(0, 0, 0, 0), "pt"),
                                                                                plot.margin = unit(c(0,0,0,0),"pt")) +
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))




#get rid of the legend for each heatmap
pval.heatmap <- pval.heatmap + theme(legend.position = "none")
mean.heatmap <- mean.heatmap + theme(legend.position = "none")

#combine each heatmap with its region label
p <- pval.heatmap + inset_element(pval.label,0,-.02,1,0,align_to = "panel",ignore_tag = TRUE)
m <- mean.heatmap + inset_element(mean.region.label,0,-.02,1,0,align_to = "panel",ignore_tag = TRUE)

#combine the two heatmaps and the barplot
bars <- wrap_elements(p + m + mean.bar) 

#to make ridgeline plots of the top 5 genera, we need to get the genus data to a long format
genus.long <- genus.rel[,colnames(genus.rel) %in% c(taxa.means$taxon[1:5],"region")] %>%
  pivot_longer(-region,names_to="taxon",values_to="rel_abundance")

#make the taxon a factor to order the faceted plots
genus.long$taxon <- factor(genus.long$taxon,levels=c(taxa.means$taxon[1:5])) 
ridges <- wrap_elements(ggplot(data=genus.long) + geom_density_ridges(aes(x=log10(rel_abundance),
                                                 y=factor(region,
                                                          levels=c('Australia/New Zealand','Central and Southern Asia',
                                                                   'Eastern and South-Eastern Asia','Latin America and the Caribbean',
                                                                   'Northern Africa and Western Asia','Sub-Saharan Africa',
                                                                   'Europe and Northern America')),
                                                 fill=region),
                                             quantile_lines = TRUE, quantiles = 2) + 
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) + theme(legend.position="none",
                                                                                axis.text.y = element_blank(),
                                                                                # axis.title.y=element_blank(),
                                                                                axis.ticks.y=element_blank() ,
                                                                                text=element_text(color='black'),
                                                                                #axis.title.x = element_blank(),
                                                                                plot.title = element_text(size=8),
                                                                                plot.margin = unit(c(10,10,10,10),"pt")
  )+ labs(x="Relative Abundance",y="Density") +  
  
  scale_x_continuous(
    labels = scales::math_format(),
    breaks = c(-5,-1)
  ) + facet_grid(~taxon)
)

#this is the layout we want for the figure
layout2 <- "
#####
AAABB
AAABB
AAABB
AAABB
CCCBB
CCCBB
CCCBB
CCCBB
"

fig4 <- overview + bars + inset_element(pval.heatmap.legend,0.15,1.02,0.59,1.07,align_to = 'plot',ignore_tag = TRUE) + 
  inset_element(mean.heatmap.legend,0.6,1.02,1,1.07,align_to = 'plot',ignore_tag = TRUE)+ ridges + plot_layout(design=layout2) +
  plot_annotation(tag_levels = 'A')
