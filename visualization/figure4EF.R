library(ggraph)
library(tidygraph)
library(igraph)

phylum_scale <- setNames(c('#BCD2EE','#832161','#06D6A0','#E88873','#6153CC','gray'),
                         c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota","other") )

edge_color_scale <- setNames(c('#8B2635','#BBDEF0'),
                             c("TRUE","FALSE") )
linetype_scale <- setNames(c('solid','dashed'),
                           c("in_ENA","not") )

#make a dataframe that contains the region names and their corresponding abbreviations
abbrs <- c("ANZ","CSA","ESEA","ENA","LAC","NAWA","SSA")
regions <- c("Australia/New Zealand","Central and Southern Asia",
             "Eastern and South-Eastern Asia","Europe and Northern America",
             "Latin America and the Caribbean", "Northern Africa and Western Asia",
             "Sub-Saharan Africa")
regionDF <<- data.frame(cbind(abbrs,regions))


#read in the p values and correlation strengths for the networks
#we'll make two networks: one for the samples within Europe/Northern America
#and the other for samples from all other regions
enaPval <- read.table('network_results/ena_pvalues.tsv',sep='\t',header=F,row.names=1)
enaCorr <- read.table('network_results/ena_corr.tsv',sep='\t',header=F,row.names=1)

otherPval <- read.table('network_results/other_pvalues.tsv',sep='\t',header=F,row.names=1)
otherCorr <- read.table('network_results/other_corr.tsv',sep='\t',header=F,row.names=1)

#we consider p > 0.001 nonsignificant
enaCorr[enaPval > 0.001] <- NA
otherCorr[otherPval > 0.001] <- NA

enaCorr[upper.tri(enaCorr,diag=T)] <- NA
otherCorr[upper.tri(otherCorr,diag=T)] <- NA

other.save <- otherCorr
ena.save <- enaCorr

#we're only considering correlations that exist in one group and NOT the other
#so if any correlation exists in the other group, set it to NA
otherCorr[!(is.na(ena.save))] <- NA
enaCorr[!(is.na(other.save))] <- NA

colnames(otherCorr) <- rownames(otherCorr)
colnames(enaCorr) <- rownames(enaCorr)

otherCorr <- otherCorr %>% rownames_to_column(var='taxon1') %>% pivot_longer(-taxon1, names_to = 'taxon2',values_to='corr', values_drop_na = TRUE) %>%
  mutate(positive = corr > 0) %>% mutate(corr = abs(corr))
ena.save <- enaCorr
colnames(enaCorr) <- rownames(enaCorr)
keep <- rownames(enaCorr) %in% append(otherCorr$taxon1,otherCorr$taxon2)
enaCorr <- enaCorr[keep,keep]
enaCorr <- enaCorr %>% rownames_to_column(var='taxon1') %>% pivot_longer(-taxon1,names_to = 'taxon2',values_to='corr',values_drop_na = FALSE) %>%
  mutate(positive = corr > 0) %>% mutate(corr = abs(corr))

otherGraph <- as_tbl_graph(otherCorr,directed = FALSE)
enaGraph <- as_tbl_graph(enaCorr,directed=F)

cpd <- readRDS(file='filtered.rds')
cpd <- make_rel(cpd)

taxa <- unique(append(otherCorr$taxon1,otherCorr$taxon2))
cpd <- cpd[,colnames(cpd) %in% taxa] #only keep the taxa we're interested in 

meta <- readRDS('metadata_from_rpackage.rds')
meta <- meta[match(rownames(cpd),meta$sample),] 
meta <- meta[meta$region %in% regionDF$regions,]
cpd <- cpd[rownames(cpd) %in% meta$sample,]

if (!(identical(rownames(cpd),meta$sample))) {
  stop('rows do not match')
}

other.cpd <- cpd[!(meta$region == "Europe and Northern America"),] 
other.cpd <- colMeans(other.cpd)
other.cpd <- other.cpd[match(names(V(otherGraph)),names(other.cpd))]
V(otherGraph)$rel_abundance <- other.cpd

ena.cpd <- cpd[meta$region == "Europe and Northern America",]
ena.cpd <- colMeans(ena.cpd)
ena.cpd <- ena.cpd[match(names(V(enaGraph)),names(ena.cpd))]
V(enaGraph)$rel_abundance <- ena.cpd
temp2 <- permute(enaGraph,match(V(enaGraph)$name,V(otherGraph)$name)) #make the order of the ENA graph match the 'other' graph

#update the taxon names for presentation
names <- read.csv('taxon_names.tsv')
names <- names[match(str_replace_all(V(otherGraph)$name,"[^[:alnum:]]",""),str_replace_all(rownames(names),"[^[:alnum:]]","")),]
if (!identical(V(otherGraph)$name,V(temp2)$name)) {
  stop('names do not match')
}
V(otherGraph)$name <- names$genus
V(temp2)$name <- names$genus

#used to scale the node sizes
max_rel_abundance <- max(append(V(otherGraph)$rel_abundance,V(temp2)$rel_abundance))
min_rel_abundance <- min(append(V(otherGraph)$rel_abundance,V(temp2)$rel_abundance))

#used to scale the edge widths
max_corr <- max(append(E(otherGraph)$corr,E(temp2)$corr),na.rm=T)
min_corr <- min(append(E(otherGraph)$corr,E(temp2)$corr),na.rm=T)

#there's not room to label all nodes, so we'll label every other node on each graph
#each node gets labeled, but only on one graph (opposite nodes are labeled on the graphs)
otherNames <- V(otherGraph)$name
otherNames[seq(from=1,to=length(otherNames),by=2)] <- ""
V(otherGraph)$label <- otherNames

#create the graph
other <- ggraph(otherGraph,layout='circle') + geom_edge_link(aes(width=corr,edge_color=positive)) + 
  geom_node_point(aes(size=rel_abundance))  + scale_color_identity(aesthetics=c('fill','color')) + labs(caption="Non-Europe and Northern America") +
  scale_color_manual(values=edge_color_scale,aesthetics=c('edge_color'),labels=c("Negative","Positive"),guide = guide_legend(title = NULL)) + 
  scale_radius(range = c(2, 5), limits = c(.001, max_rel_abundance),breaks=c(0.025,0.075,0.125),name="Relative Abundance",
               guide=guide_legend(title.position="bottom")) + 
  scale_edge_width(range=c(1,5),limits = c(0.01,max_corr),breaks=c(0.05,0.15),name="Correlation",guide=guide_legend(title.position="bottom")) +
  theme(legend.position="bottom", plot.caption = element_text(hjust=0.5,size=10))

#add another layer of the data so we can get rid of the pieces we don't need 
#we will use this to label the nodes with geom_segment
other$data2 <- other$data
other$data2[seq(from=1,to=length(otherNames),by=2),] <- NA

#add the taxa names
other <- other + 
  scale_x_continuous(expand=expansion(mult=0.2)) + scale_y_continuous(expand=expansion(mult=0.2)) + 
  geom_segment(aes(x=other$data2$x*1.05,xend=other$data2$x*1.25,y=other$data2$y*1.05,yend=other$data2$y*1.25)) + 
  geom_node_label(aes(label=stringr::str_wrap(label,10),lineheight=0.7),nudge_x = other$data$x * .75, nudge_y = other$data$y * .3,size=2.5,fill='white',color='white') +
  geom_node_text(aes(label=stringr::str_wrap(label,10),lineheight=0.7),nudge_x = other$data$x * .75, nudge_y = other$data$y * .3,size=2.5)

network.leg <- ggpubr::as_ggplot(ggpubr::get_legend(other))

other <- other + theme(legend.position="none")

enaNames <- V(temp2)$name
enaNames[seq(from=2,to=length(otherNames),by=2)] <- ""
V(temp2)$label <- enaNames

ena <- ggraph(temp2,layout='circle') + geom_edge_link(aes(width=corr,edge_color=positive)) + 
  geom_node_point(aes(size=rel_abundance))  + scale_color_identity(aesthetics=c('fill','color')) + labs(caption="Europe and Northern America") +
  scale_color_manual(values=edge_color_scale,aesthetics=c('edge_color')) + scale_radius(range = c(2, 5), limits = c(.001, max_rel_abundance)) + 
  scale_edge_width(range=c(1,5),limits = c(0.01,max_corr)) + theme(legend.position="none", plot.caption = element_text(hjust=0.5,size=10))

ena$data2 <- ena$data
ena$data2[seq(from=2,to=length(otherNames),by=2),] <- NA

ena <- ena  + 
  scale_x_continuous(expand=expansion(mult=0.2)) + scale_y_continuous(expand=expansion(mult=0.2)) + 
  geom_segment(aes(x=ena$data2$x*1.05,xend=ena$data2$x*1.25,y=ena$data2$y*1.05,yend=ena$data2$y*1.25)) + 
geom_node_label(aes(label=stringr::str_wrap(label,10),lineheight=0.7),nudge_x = ena$data$x * .75, nudge_y = ena$data$y * .3,size=2.5,fill='white',color='white') +
geom_node_text(aes(label=stringr::str_wrap(label,10),lineheight=0.7),nudge_x = ena$data$x * .75, nudge_y = ena$data$y * .3,size=2.5)

networks <- wrap_plots(ena + other)
