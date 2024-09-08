library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(tidyr) # pivot_longer
library(magrittr) # for and()
library(vegan)
library(grid) # textGrob

#########
#Figure 1C, Figure 1D, Figure 1E
#########
prevalence.phylum <- get_prevalence(readRDS('phylum.rds'))
prevalence.class <- get_prevalence(readRDS('class.rds'))
prevalence.order <- get_prevalence(readRDS('order.rds'))
prevalence.family <- get_prevalence(readRDS('family.rds'))
prevalence.genus <- get_prevalence(readRDS('filtered.rds'))
#---------- Most prevalent taxa

taxa_per_chart <- 10

# We leave the top 5 most common phyla, and bunch
# the rest together in an "other" bucket
top_phyla <- prevalence.phylum
top_phyla$taxon <- gsub('^\\w+\\.(\\w+)$', '\\1', top_phyla$taxon)
top_phyla <- top_phyla %>%
  arrange(desc(samples)) %>%
  slice_head(n=5) %>%
  rbind(data.frame(
    taxon=c('other'),
    samples=c(0)
  ))


### Phylum level
toplot.phylum <- prevalence.phylum %>%
  arrange(desc(samples)) %>%
  slice_head(n=taxa_per_chart)

# Standardize taxa names
toplot.phylum$phylum <- gsub('^\\w+\\.(\\w+)$', '\\1', toplot.phylum$taxon)
# only keep names for the top 5 phyla
toplot.phylum$phylum <- if_else(toplot.phylum$phylum %in% top_phyla$taxon, toplot.phylum$phylum, 'other')

toplot.phylum$taxon <- format_names(toplot.phylum$taxon)

### Class level
toplot.class <- prevalence.class %>%
  arrange(desc(samples)) %>%
  slice_head(n=taxa_per_chart)

toplot.class$phylum <- gsub('^\\w+\\.(\\w+)\\..+$', '\\1', toplot.class$taxon)
toplot.class$phylum <- if_else(toplot.class$phylum %in% top_phyla$taxon, toplot.class$phylum, 'other')
toplot.class$taxon <- format_names(toplot.class$taxon)

### Order level
toplot.order <- prevalence.order %>%
  arrange(desc(samples)) %>%
  slice_head(n=taxa_per_chart)

# Annotate each order with its phylum
toplot.order$phylum <- gsub('^\\w+\\.(\\w+)\\..+$', '\\1', toplot.order$taxon)
toplot.order$phylum <- if_else(toplot.order$phylum %in% top_phyla$taxon, toplot.order$phylum, 'other')

toplot.order$taxon <- format_names(toplot.order$taxon)

### Family level
toplot.family <- prevalence.family %>%
  arrange(desc(samples)) %>%
  slice_head(n=taxa_per_chart)

toplot.family$phylum <- gsub('^\\w+\\.(\\w+)\\..+$', '\\1', toplot.family$taxon)
toplot.family$phylum <- if_else(toplot.family$phylum %in% top_phyla$taxon, toplot.family$phylum, 'other')
toplot.family$taxon <- format_names(toplot.family$taxon)

### Genus level
toplot.genus <- prevalence.genus %>%
  arrange(desc(samples)) %>%
  slice_head(n=taxa_per_chart)
toplot.genus$phylum <- gsub('^\\w+\\.(\\w+)\\..+$', '\\1', toplot.genus$taxon)
toplot.genus$phylum <- if_else(toplot.genus$phylum %in% top_phyla$taxon, toplot.genus$phylum, 'other')
toplot.genus$taxon <- format_names(toplot.genus$taxon)
toplot.genus$taxon[10] <- '[Ruminococcus] torques group' # fixing funny formatting
######## PLOTTING
# Color scheme
# https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7

# Convert to latest phylum names for paper
update_phyla <- function(df, change=TRUE) {
  df[df$phylum=='Firmicutes',]$phylum <- 'Bacillota'
  df[df$phylum=='Proteobacteria',]$phylum <- 'Pseudomonadota'
  df[df$phylum=='Actinobacteriota',]$phylum <- 'Actinomycetota'
  if(change) df$taxon <- gsub('^.+ (\\w+)$', '\\1', df$taxon)
  return(df)
}

toplot.phylum <- update_phyla(toplot.phylum)
toplot.class<- update_phyla(toplot.class)
toplot.order <- update_phyla(toplot.order)
toplot.family <- update_phyla(toplot.family)
# we change the taxon names differently for genera because some
# are multiple words:
toplot.genus <- update_phyla(toplot.genus, change=FALSE)
toplot.genus$taxon <- gsub('^(?:[^ ]+ ){4}(.+)$', '\\1', toplot.genus$taxon)


phylum_palette <- setNames(c('#BCD2EE','#832161','#06D6A0','#E88873','#6153CC','gray'),
                         c("Bacillota","Pseudomonadota","Actinomycetota","Bacteroidota", "Desulfobacterota","other") )
phylum_scale <- scale_fill_manual(values=phylum_palette, aesthetics=c('colour','fill'))

# shared axis:
number_scale <- scale_y_continuous(
  labels = label_number(scale_cut = cut_short_scale()),
  limits = c(0,155000),
  breaks = (c(0,50000,100000,150000)),
  expand=c(0,0)
)

othersize = 3.8 # element_text is 11, what's the number for everything else?

label_geom <- geom_text(aes(x=taxon, y=1000, label=taxon, color=NA),
                        hjust=0, size=othersize, color='black')

shared_theme <- theme(
  legend.position='none',
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  axis.text.y = element_blank(),
  axis.text.x = element_text(size=11),
  axis.ticks.y = element_blank(),
  plot.margin = margin(0, 4, 16, 4, "pt"),
  rect = element_rect(fill='transparent')
)

# Adjust outdated names
toplot.phylum[toplot.phylum$taxon=='Firmicutes',]$taxon <- 'Bacillota'
toplot.phylum[toplot.phylum$taxon=='Proteobacteria',]$taxon <- 'Pseudomonadota'
toplot.phylum[toplot.phylum$taxon=='Actinobacteriota',]$taxon <- 'Actinomycetota'

plot.phylum <- ggplot(toplot.phylum, aes(x=reorder(taxon, samples),
                                         y=samples, fill=phylum, color=phylum)) +
  geom_bar(stat="identity", alpha=0.5) +
  label_geom +
  theme_bw() +
  coord_flip() +
  number_scale +
  labs(x='Phylum', y='') +
  phylum_scale +
  shared_theme +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=11, margin = margin(r = -2.1, unit = "cm"))
  )

tagadjust = c(0.05, 0.965) # scooting labels in the middle of the second row

plot.class <- ggplot(toplot.class, aes(x=reorder(taxon, samples),
                                       y=samples, fill=phylum, color=phylum)) +
  geom_bar(stat="identity", alpha=0.5) +
  label_geom +
  theme_bw() +
  number_scale +
  coord_flip() +
  labs(x='Class', y='Samples') +
  phylum_scale +
  shared_theme +
  theme(
    plot.tag.position = tagadjust
  )

plot.order <- ggplot(toplot.order, aes(x=reorder(taxon, samples),
                                       y=samples, fill=phylum, color=phylum)) +
  geom_bar(stat="identity", alpha=0.5) +
  label_geom +
  theme_bw() +
  number_scale +
  coord_flip() +
  labs(x='Order', y='') +
  phylum_scale +
  shared_theme +
  theme(
    plot.tag.position = tagadjust
  )

##### Supplementary panels at lower taxonomic levels
plot.family <- ggplot(toplot.family, aes(x=reorder(taxon, samples),
                                       y=samples, fill=phylum, color=phylum)) +
  geom_bar(stat="identity", alpha=0.5) +
  label_geom +
  theme_bw() +
  number_scale +
  coord_flip() +
  labs(x='Family', y='', fill='Phylum', color='Phylum') +
  phylum_scale +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=11),
    axis.ticks.y = element_blank()
  )

plot.genus <- ggplot(toplot.genus, aes(x=reorder(taxon, samples),
                                       y=samples, fill=phylum, color=phylum)) +
  geom_bar(stat="identity", alpha=0.5) +
  label_geom +
  theme_bw() +
  number_scale +
  coord_flip() +
  labs(x='Genus', y='Prevalence (samples)') +
  phylum_scale +
  theme(
    legend.position='none',
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=11),
    axis.ticks.y = element_blank()
  )

plot.phylum / plot.class / plot.order

plot.family / plot.genus

ggsave('figures/panels/supp_prevalence.pdf', device='pdf', height=9, width=9, units='in')

# Supplementary table
write.table(rbind(prevalence.phylum, prevalence.class, prevalence.order, prevalence.family, prevalence.genus),
          file='figures/panels/prevalence.tsv', row.names=F, sep='\t')

##########################################
# RAREFACTION
##########################################
# These results are generated by rarefaction_cluster.R
rarefaction <- readRDS('rarefaction.rds') %>%
  group_by(level, scount) %>%
  summarise(
    sd = sd(observed, na.rm = TRUE),
    observed = mean(observed)
  )

rarefaction$level <- as.character(rarefaction$level)
rarefaction[rarefaction$level=='genus',]$level <- 'Genus'
rarefaction[rarefaction$level=='family',]$level <- 'Family'
rarefaction[rarefaction$level=='order',]$level <- 'Order'
rarefaction[rarefaction$level=='class',]$level <- 'Class'
rarefaction[rarefaction$level=='phylum',]$level <- 'Phylum'
rarefaction$level <- as.factor(rarefaction$level)

rarefaction.labels <- rarefaction %>%
  filter(scount == max(rarefaction$scount)) %>%
  mutate(scount=100000, observed=observed+200)

# manual fix for genus
rarefaction.labels[rarefaction.labels$level=='Genus',]$observed <- 3500
rarefaction.labels[rarefaction.labels$level=='Phylum',]$observed <- -80

rare_plot <- ggplot(rarefaction, aes(x=scount, y=observed,
                    ymin = observed-sd, ymax = observed+sd,
                    color=level, label=level)) +
  geom_errorbar(width = 3000) +
  geom_text(data=rarefaction.labels, aes(x=scount, y=observed, color=level, label=level),
            fontface='bold') +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_brewer(palette='Dark2') +
  scale_x_continuous(
    breaks=c(0,50000,100000,150000),
    labels=label_number(scale_cut = cut_short_scale())
  ) +
  scale_y_continuous(labels=scales::comma, expand=expansion(add=400,0)) +
  theme_bw() +
  labs(x='Samples', y='Unique taxa', color='Taxonomic level') +
  theme(
    legend.position = 'none',
    axis.text = element_text(size=11)
  )

###### RAREFACTION ANALYSIS
maxes <- rarefaction %>%
  group_by(level) %>%
  summarise(max=max(observed))

# calculate discovery rate
rate <- rarefaction %>%
  mutate(prev = ifelse(scount > 1, lag(observed), 0)) %>%
  mutate(rate = (observed-prev)/scount) %>%
  inner_join(maxes, by='level') %>%
  mutate(rate_of_max = rate / max) %>%
  mutate(sample_per = 1/rate)
################################

## PHYLUM HISTOGRAM
taxphylum <- readRDS('phylum.rds')

final.rel <- make_rel(taxphylum)
colnames(final.rel) <- gsub('^Bacteria\\.([^\\(])', '\\1', colnames(final.rel))

# make histograms
full.long <- final.rel %>%
  pivot_longer(everything(), names_to = "phylum", values_to = "rel") %>%
  filter(phylum %in% top_phyla$taxon) %>%
  update_phyla(change=FALSE)

# Hard-coding phylum names here now that we know what the top ones are
top_levels <- c("Bacillota", "Pseudomonadota", "Actinomycetota", "Bacteroidota", "Desulfobacterota", "other")

full.long$phylum <- factor(full.long$phylum, levels=top_levels)


rel_histogram <- ggplot(full.long, aes(x=rel, y=after_stat(count), color=phylum)) +
  geom_freqpoly(linewidth=1) +
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks=c(100, 1000, 1e4, 1e5)
  ) +
  scale_x_continuous(labels = scales::percent, expand=c(0,0)) +
  coord_cartesian(xlim=c(0,1)) +
  phylum_scale +
  theme_bw() +
  labs(
    x='Relative abundance',
    y='Samples',
    color='Taxon'
  ) +
  theme(
    axis.text = element_text(size=11),
    legend.text = element_text(size=11),
    legend.position = 'none'
  )

linear <- ggplot(full.long, aes(x=rel, y=after_stat(count), color=phylum)) +
  geom_freqpoly(linewidth=1) +
  scale_y_continuous(
    #labels = scales::label_number(scale_cut = scales::cut_short_scale()), # doesn't work in 2024 for some reaon
    expand=c(0,0)
  ) +
  scale_x_continuous(labels = scales::percent, expand=c(0,0)) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,10000)) +
  theme_bw() +
  phylum_scale +
  labs(
    x='Relative abundance',
    y='Samples',
    color='Taxon'
  ) +
  theme(
    axis.text = element_text(size=11),
    legend.text = element_text(size=11)
  )
ggsave(plot=linear, 'figures/panels/supp_relhistogram.pdf', device='pdf', height=6, width=9, units='in')


## DIVERSITY HISTOGRAM
divdata <- readRDS('filtered.rds')

shannon <- diversity(divdata, "shannon")
shannon <- as.data.frame(shannon)
shannon$srr <- rownames(shannon)

depth <- data.frame(rownames(divdata), rowSums(divdata))
colnames(depth) <- c('srr','depth')

diversity <- ggplot(shannon, aes(x=shannon)) +
  geom_histogram(bins=50, fill='#000000', color='#000000') +
  geom_vline(xintercept=median(shannon$shannon),
             linewidth=1,color="gray") +
  annotate("text", label=paste(
      "median:", round(median(shannon$shannon), 2)
    ),
    size = othersize,
    x=2.3,y=8000#x=3.5, y=9000 asdf
  ) +
  theme_bw() +
  labs(
    x="Shannon diversity",
    y="Samples"
  ) +
  scale_y_continuous(labels=scales::comma) +
  theme(
    axis.text = element_text(size=11),
    rect = element_rect(fill='transparent')
  )

# characterizing the distribution
nrow(depth)
dtest <- depth %>% filter(depth >= 10000) %>%
  filter(depth <= 45000)

(nrow(dtest) / nrow(depth))*100

dtest <- depth %>% filter(depth > 85000)
(nrow(dtest) / nrow(depth))*100

depth_plot <- ggplot(depth, aes(x=depth)) +
  geom_histogram(bins=50, fill='#000000', color='#000000') +
  geom_vline(xintercept=median(depth$depth),
             linewidth=1,color="gray") +
  annotate("text", label=paste(
      "median:", comma(median(depth$depth))
    ),
    size = othersize,
    x=7e5,y=9000#x=500000, y=9000
  ) +
  theme_bw() +
  labs(
    x="Merged reads",
    y="Samples"
  ) +
  scale_y_continuous(
    labels=scales::comma,
    limits=c(0,11250),
    expand=c(0,0)
  ) +
  scale_x_log10(labels=scales::label_log()) +
  theme(
    axis.text=element_text(size=11)
  )


###########################
#------- assembling the whole figure

row1 <- (wrap_elements(panel = textGrob('Pipeline diagram here')) | depth_plot) +
  plot_layout(widths=c(2,1))

fullfig <- row1 /
  (plot.phylum + plot.class + plot.order) /
  (rel_histogram + diversity + rare_plot ) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag=element_text(face='bold')
  )

ggsave('figures/panels/figure1.pdf', plot=fullfig,
        device='pdf', height=10, width=10, units='in')
