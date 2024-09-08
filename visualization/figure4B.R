#color scheme used for the regions
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

#color scheme used for the phyla
phylum_scale <- setNames(c('#BCD2EE','#832161','#06D6A0','#E88873','#6153CC','gray'),
                         c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota","other") )


res <- readRDS("subsample_rarefied_results.rds")

#grab the total number of taxa found in each iteration and create a separate dataframe
totalTaxa <- res[res$stat == "total_taxa",]
totalTaxa <- totalTaxa %>% select(-stat)  %>% pivot_longer(-region,names_to="iter",values_to = "total_taxa")
  

res <- res %>% pivot_longer(cols=-c(region,stat),names_to="iter",values_to="mean")

res <- res %>% full_join(totalTaxa,by=c("region","iter"))

#remove the 'total_taxa' stat
remove <- res$stat == "total_taxa"
res <- res[-remove,]
res <- res %>% separate_wider_delim(stat,delim="_",names=c('stat','level'))

#grab only the info for NA taxa (not unID taxa)
meanTaxa <- res[res$stat == "naTaxa",]

#add a stat for the proportion of NA taxa
meanTaxa$prop <- meanTaxa$mean / meanTaxa$total_taxa

#calculate the mean and standard deviation for the mean number of taxa discovered for each group (region, taxonomic rank)
meanTaxa <- meanTaxa %>% group_by(region,level) %>% summarise(meanTaxa = mean(mean),sdTaxa = sd(mean),meanProp = mean(prop),sdProp=sd(prop))


colError1 <- ggplot(meanTaxa) + geom_col(aes(x=factor(level,levels=c("genus","family","order","class","phylum")),y=meanTaxa,fill=region,group=region),position="dodge")  + 
  geom_errorbar(aes(x=factor(level,levels=c("genus","family","order","class","phylum")),ymin=meanTaxa-sdTaxa,ymax=meanTaxa+sdTaxa,group=region),position = "dodge") + 
  scale_fill_manual(values=region_scale, aesthetics=c('colour','fill')) + 
  ylab('Average number of unidentified taxa')  + theme_bw() +theme(axis.title.x=element_blank(), legend.position = 'none',
                                                                   axis.text.x = element_text(angle=45,hjust=1)
  ) 
