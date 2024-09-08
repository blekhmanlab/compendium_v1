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
meta <- readRDS('metadata_from_rpackage.rds')

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
  theme(plot.margin = unit(c(5,15,0,0),"pt"),
        axis.title = element_text(size=11),
        legend.position="bottom",
        legend.key.spacing.y = unit(1, "pt"),
        legend.box.spacing = unit(1,"pt")) +
  labs(
    x='Relative abundance',
    y='Samples'
  ) + 
  facet_wrap(~region,nrow=1, labeller = labeller(region = label_wrap_gen(13))) + 
  scale_fill_manual(values=phylum_scale, aesthetics=c('colour','fill')) + 
  guides(color=guide_legend(title="Phylum"))
