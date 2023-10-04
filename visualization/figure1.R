library(ggplot2)
library(patchwork)
library(scales)
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

taxa_per_chart <- 10 # how many in each panel?

# For the color scheme, we leave the top 5 most common phyla
# and bunch the rest together in an "other" bucket
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
# Make display names pretty
toplot.phylum$taxon <- format_names(toplot.phylum$taxon)

######## Convert new names back to the old ones:
toplot.phylum[toplot.phylum$taxon=='Actinobacteriota',]$taxon <- 'Actinobacteria'
toplot.phylum[toplot.phylum$taxon=='Bacteroidota',]$taxon <- 'Bacteroidetes'
toplot.phylum[toplot.phylum$taxon=='Verrucomicrobiota',]$taxon <- 'Verrucomicrobia'
toplot.phylum[toplot.phylum$taxon=='Fusobacteriota',]$taxon <- 'Fusobacteria'


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

######## PLOTTING
# Color scheme
# https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7

phylum_palette <- setNames(c('#BCD2EE','#832161','#06D6A0','#E88873','#6153CC','gray'),
                         c("Firmicutes","Proteobacteria","Actinobacteriota","Bacteroidota", "Desulfobacterota","other") )
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

plot.phylum <- ggplot(toplot.phylum, aes(x=reorder(taxon, samples),
                                         y=samples, fill=phylum, color=phylum)) +
  geom_bar(stat="identity", alpha=0.5) +
  #geom_text(aes(x=taxon, y=1000, label=taxon), hjust=0, size=4, color='black') +
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

suppfig <- plot.family / plot.genus

ggsave('supp_prevalence.pdf', plot=suppfig, device='pdf',
       height=9, width=9, units='in'
)

# Supplementary table 1
write.table(rbind(
    prevalence.phylum,
    prevalence.class,
    prevalence.order,
    prevalence.family,
    prevalence.genus
  ),
  file='prevalence.tsv', row.names=F, sep='\t'
)

##########################################
# Figure 1F, Figure 1G: RELATIVE ABUNDANCE
##########################################
taxphylum <- readRDS('phylum.rds')

# Start with the phylum-level data, and turn
# it into a relative abundance table
final.rel <- make_rel(taxphylum)
colnames(final.rel) <- gsub('^Bacteria\\.([^\\(])', '\\1', colnames(final.rel))
# we can't possibly display 150,000+ samples on a screen
# that's MAYBE 3000 pixels wide, so grab a random sample
# of 5000:
set.seed(43)
subset <- final.rel[sample(nrow(final.rel), 5000), ]

# if we want to order the samples using certain taxa,
# we need an extra copy of those, which we'll attach
# to every entry for every sample in the pivot_longer version.
secondmax.name <- function(x) {
  names <- data.frame(cols=colnames(subset)) %>%
    filter(!cols %in% names(which.max(x)))

  names$cols[which.max(x[x != max(x)])]
}
secondmax.val <- function(x) {
  max(x[x != max(x)])
}

# annotate each sample with most prevalent taxon
plurality_values <- colnames(subset)[apply(subset,1,which.max)]
# annotate each sample with SECOND-MOST prevalent taxon
second_values <- apply(subset,1,secondmax.name)

subset$sorter1 <- apply(subset,1,max)
subset$sorter2 <- apply(subset,1,secondmax.val)
subset$plurality1 <- plurality_values
subset$plurality2 <- second_values

subset$sample <- rownames(subset)
# Only keep the sorting values if they're in the top 5 most
# prevalent phyla
subset[!(subset$plurality1 %in% top_phyla$taxon),]$plurality1 <- 'other'
subset[!(subset$plurality2 %in% top_phyla$taxon),]$plurality2 <- 'other'

# now go to long form to plot
final.long <- subset %>%
  pivot_longer(!c(sample, plurality1, plurality2, sorter1,sorter2),
               names_to = "taxon", values_to = "rel")

final.long[!(final.long$taxon %in% top_phyla$taxon),]$taxon <- 'other'

# Sort the samples using the abundance of the most prevalent
# taxon, followed by the total abundance of the TWO most prevalent
# taxa
final.long$sample <- factor(final.long$sample, levels=subset[
  rev(order(
    subset$sorter1,
    subset$sorter1+subset$sorter2
  ))
  ,]$sample)
final.long$plurality1 <- factor(final.long$plurality1, levels=top_phyla$taxon)
final.long$plurality2 <- factor(final.long$plurality2, levels=top_phyla$taxon)
final.long$taxon <- factor(final.long$taxon, levels=top_phyla$taxon)

stacked <- ggplot(final.long, aes(fill=taxon, y=rel, x=sample)) +
  geom_bar(stat="identity",width=1) +
  phylum_scale +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  labs(x='Sample',y='Relative abundance', fill='Phylum') +
  theme_bw() +
  theme(
    legend.position='none',
    axis.text.y = element_text(size=11),
    axis.title.y = element_text(size=11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  facet_grid(
    cols=vars(plurality1,plurality2),
    space='free_x',
    scales='free_x',
    drop=T
  )

# Supplementary Figure 2
# (conventionally ordered stacked bar)
final.long$sample <- factor(final.long$sample, levels=subset[
  order(subset[[top_phyla$taxon[1]]],
        subset[[top_phyla$taxon[1]]]+subset[[top_phyla$taxon[2]]],
        subset[[top_phyla$taxon[1]]]+subset[[top_phyla$taxon[2]]]+subset[[top_phyla$taxon[3]]]
    ),
]$sample)

stacked.conventional <- ggplot(final.long, aes(fill=taxon, y=rel, x=sample)) +
  geom_bar(stat="identity",width=1) +
  phylum_scale +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  labs(x='Sample',y='Relative abundance', fill='Phylum') +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=11),
    axis.title.y = element_text(size=11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(plot=stacked.conventional, 'supp_stacked.pdf', device='pdf', height=3.6, width=10, units='in')

### Supplementary table 2
#   Evaluate patterns for the filtered dataset
combos <- final.rel
full.plurality_values <- colnames(combos)[apply(combos,1,which.max)]
full.second_values <- apply(combos,1,secondmax.name)

combos$sorter1 <- apply(combos,1,max)
combos$sorter2 <- apply(combos,1,secondmax.val)
combos$plurality1 <- full.plurality_values
combos$plurality2 <- full.second_values

combos$sample <- rownames(combos)
combos[!(combos$plurality1 %in% top_phyla$taxon),]$plurality1 <- 'other'
combos[!(combos$plurality2 %in% top_phyla$taxon),]$plurality2 <- 'other'

# determining frequency of combinations
combos$toptwo <- combos$sorter1 + combos$sorter2

# how many have firmicutes in the top two?
combos %>%
  filter(or(
      plurality1 == 'Firmicutes',
      plurality2 == 'Firmicutes'
    )
  ) %>%
  nrow()

# Supplementary Table 2
write.table(table(combos$plurality1, combos$plurality2), file='panels/combos.tsv', sep='\t')

#########
# Figure 1G
#########

# make histograms
full.long <- final.rel %>%
  pivot_longer(everything(), names_to = "taxon", values_to = "rel") %>%
  filter(taxon %in% top_phyla$taxon)

full.long$taxon <- factor(full.long$taxon, levels=c(top_phyla$taxon))


rel_histogram <- ggplot(full.long, aes(x=rel, y=after_stat(count), color=taxon)) +
  geom_freqpoly(linewidth=1) +
  #scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks=c(100, 1000, 1e4, 1e5)
  ) +
  scale_x_continuous(labels = scales::percent, expand=c(0,0)) +
  coord_cartesian(xlim=c(0,1)) +
  #scale_fill_manual(values=old_names, aesthetics=c('colour','fill')) +
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

# Supplementary Figure 3
linear <- ggplot(full.long, aes(x=rel, y=after_stat(count), color=taxon)) +
  geom_freqpoly(linewidth=1) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    expand=c(0,0)
  ) +
  scale_x_continuous(labels = scales::percent, expand=c(0,0)) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,10000)) +
  #scale_fill_manual(values=old_names, aesthetics=c('colour','fill')) +
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
ggsave(plot=linear, 'supp_relhistogram.pdf', device='pdf', height=6, width=9, units='in')

#########
# Figure 1H
#########

divdata <- readRDS('filtered.rds')

shannon <- diversity(divdata, "shannon") %>%
  as.data.frame() %>%
  dplyr::rename(shannon='.')
shannon$srr <- rownames(shannon)

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

#########
# Figure 1B
#########

depth <- data.frame(rownames(divdata), rowSums(divdata))
colnames(depth) <- c('srr','depth')

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


##########################################
# Figure 1I, RAREFACTION
##########################################

rarefaction <- readRDS('rarefaction.rds') %>%
  group_by(level, scount) %>%
  summarise(
    sd = sd(observed, na.rm = TRUE),
    observed = mean(observed)
  )

# Fix the labels
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

# manual fix for labels
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
  #scale_x_log10(labels = comma_format()) +
  theme(
    #legend.text = element_text(size=NA),
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

#------- assembling the whole figure

row1 <- (wrap_elements(panel = textGrob('Pipeline diagram here')) | depth_plot) +
  plot_layout(widths=c(2,1))

fullfig <- row1 /
  (plot.phylum + plot.class + plot.order) /
  stacked /
  (rel_histogram + diversity + rare_plot ) +
plot_annotation(tag_levels = 'A') &
theme(
  plot.tag=element_text(face='bold')
  #plot.margin = margin(0, 10, 0, 5, "pt")
)

ggsave('panels/figure1.pdf', plot=fullfig,
       device='pdf', height=10.5, width=8, units='in')
