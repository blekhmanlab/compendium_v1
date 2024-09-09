library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(tidyr) # for replace_na
library(ggnewscale)

figtheme <- theme(
  axis.text = element_text(size=11),
  axis.title = element_text(size=12)
)

### ASSAYS OVER TIME
###############################
libs = c('AMPLICON','WGS','OTHER','WGA','RNA-Seq')

libcolors <- setNames(c(
  "#ca0020",
  "#8ad1e9",
  "#f4a582",
  "#929191",
  "#404040"
), libs)
# colors:
libscale <- scale_fill_manual(values=libcolors,
                              aesthetics=c('fill', 'color'))

data <- read.delim('sra_samples.tsv')

data$year <- strtoi(substr(data$pubdate, 1, 4))

data <- data  %>%
  filter(project != 'PRJNA803937') %>% # single-cell paper w/ 50k samples
  filter(year < 2024) %>%
  filter(year > 2010)

find_biggest <- data %>%
  group_by(library_strategy) %>%
  summarise(samples=n()) %>%
  top_n(5, samples)

supptable <- data %>%
  group_by(library_strategy, year) %>%
  summarise(
    samples=n(),
    projects=n_distinct(project)
  )
write.csv(supptable, file='figures/assaycounts.csv',row.names=F)

yearcatcounts <- data %>%
  mutate(keep = library_strategy %in% find_biggest$library_strategy) %>%
  mutate(plotcat = ifelse(keep, library_strategy, 'OTHER')) %>%
  group_by(plotcat, year) %>%
  summarise(
    samples=n(),
    projects=n_distinct(project)
  ) %>%
  mutate(plotcat = factor(plotcat,
                          levels=libs)
  ) %>%
  filter(year < 2024) %>%
  filter(year >= 2010) %>%
  mutate(running_samples = cumsum(samples)) %>%
  mutate(running_projects = cumsum(projects))

assaytime <- ggplot(yearcatcounts, aes(x=year, y=running_samples, color=plotcat)) +
  geom_line(linewidth=1) +
  scale_y_continuous(labels=scales::comma) +
  libscale +
  labs(x='Year', y='Cumulative samples', color='Library strategy') +
  theme_bw() +
  figtheme +
  theme(
    legend.title = element_text(size=9),
    legend.text = element_text(size=9)
  )

library(cowplot) # for get_legend

assaylegend <- get_legend(assaytime)

assaycombined <- assaytime + theme(legend.position='none') +
  inset_element(assaylegend, 0.34, 0.75, 0.5, 0.55,
    align_to='plot',
    on_top=T
  )


# Supplemental figure: assay counts by PROJECT
assayprojects <- ggplot(yearcatcounts, aes(x=year, y=running_projects, color=plotcat)) +
    geom_line(linewidth=1) +
    scale_y_continuous(labels=scales::comma) +
    libscale +
    labs(x='Year', y='Cumulative projects', color='Strategy') +
    theme_bw() +
    figtheme +
    theme(
      legend.position = 'none'
    )

yearcatcounts$plotcat_reverse <- factor(yearcatcounts$plotcat, levels=rev(libs))
assayprojects.year <- ggplot(yearcatcounts, aes(x=year, y=projects, fill=plotcat_reverse)) +
  geom_bar(stat='identity') +
  scale_y_continuous(labels=scales::comma) +
  scale_fill_manual( # this is the "libscale" with the order flipped back to normal
    values=libcolors,
    aesthetics=c('fill', 'color'),
    guide = guide_legend(reverse = TRUE, nrow=2)
  ) +
  labs(x='Year', y='Projects in year', fill='Strategy') +
  theme_bw() +
  figtheme +
  theme(
    legend.position = 'top',
    legend.text = element_text(size=8)
  )

assay_counts <- yearcatcounts %>%
  group_by(plotcat) %>%
  summarise(
    samples=sum(samples)
  )
## Samples per instrument

biggest_inst <- data %>%
  filter(library_strategy %in% c('WGS','AMPLICON')) %>%
  group_by(library_strategy, instrument) %>%
  summarise(samples=n()) %>%
  ungroup() %>%
  group_by(library_strategy) %>%
  slice_max(samples, n=7)

instcounts <- data %>%
  filter(library_strategy %in% c('WGS','AMPLICON')) %>%
  mutate(keep = paste0(library_strategy,instrument) %in% paste0(biggest_inst$library_strategy,biggest_inst$instrument)) %>%
  mutate(plotcat = ifelse(keep, instrument, 'OTHER')) %>%
  group_by(library_strategy, plotcat, year) %>%
  summarise(
    samples=n(),
    projects=n_distinct(project)
  ) %>%
  filter(year < 2024) %>%
  filter(year >= 2010) %>%
  mutate(running_samples = cumsum(samples)) %>%
  mutate(running_projects = cumsum(projects))

inst_top <- ggplot(instcounts[instcounts$library_strategy=='AMPLICON',], aes(x=year, y=running_samples, color=plotcat)) +
  geom_line(linewidth=1) +
  scale_y_continuous(labels=scales::comma) +
  labs(x='Year', y='Cumulative samples', color='Instrument', title='Amplicon') +
  theme_bw() +
  figtheme +
  theme(
    legend.title = element_text(size=9),
    legend.text = element_text(size=9),
    legend.position = 'right'
  )

inst_bottom <- ggplot(instcounts[instcounts$library_strategy=='WGS',], aes(x=year, y=running_samples, color=plotcat)) +
  geom_line(linewidth=1) +
  scale_y_continuous(labels=scales::comma) +
  labs(x='Year', y='Cumulative samples', color='Instrument', title='Shotgun') +
  theme_bw() +
  figtheme +
  theme(
    legend.title = element_text(size=9),
    legend.text = element_text(size=9),
    legend.position = 'right'
  )

# SUPPLEMENTARY FIGURE 15
towrite <- assayprojects + inst_top +
  assayprojects.year + inst_bottom +
  plot_annotation(tag_levels='A')

ggsave(plot=towrite,
  file='supp_tech.pdf',
  device='pdf', height=9, width=11, units='in')

tempcounting <- instcounts %>%
  filter(library_strategy=='AMPLICON', year == 2023)
sum(tempcounting$running_samples)

## Assays by region
####################

### Assays by REGION rather than time:
regionlevels <- data %>%
  filter(region!='') %>%
  group_by(region) %>%
  summarise(samples=n()) %>%
  arrange(samples)

regioncounts <- data %>%
  filter(region!='') %>%
  filter(library_strategy %in% c('AMPLICON','WGS')) %>%
  group_by(region, library_strategy) %>%
  mutate(region=factor(region, levels=regionlevels$region)) %>%
  mutate(library_strategy=factor(library_strategy, levels=rev(c('AMPLICON','WGS')))) %>%
  summarise(samples=n())

regiontotals <- regioncounts %>%
  group_by(region) %>%
  summarise(samples=sum(samples))

# make plots using only 16S and shotgun:
proportional <- ggplot(regioncounts[regioncounts$library_strategy!='OTHER',], aes(y=region, x=samples, fill=library_strategy)) +
  geom_bar(stat='identity', position='fill') +
  libscale +
  scale_x_continuous(breaks=c(0,0.5, 1)) +
  labs(x='Proportion', y='Region') +
  theme_bw() +
  figtheme +
  theme(
    legend.position = 'none'
  )

# Position of the total sample labels
labelx <-  + c(20000, rep(90000, 5), 125000, 215000, 310000)

raw <- ggplot(regioncounts, aes(x=samples, y=region)) +
  geom_bar(mapping=aes(
    x=samples,
    y=region,
    fill=library_strategy
  ), stat='identity') +
  geom_text(regiontotals, mapping=aes(
    x=labelx, y=region, label=scales::comma(samples))
  ) +
  libscale +
  scale_x_continuous(
    labels = unit_format(unit = 'K',scale = 1e-3, sep=''),
    breaks=c(1e4, 2e5, 4e5)
  ) +
  labs(
    x='Total samples', y='',
    fill='Library strategy'
  ) +
  theme_bw() +
  figtheme +
  # colors listed in wrong order for some reason
  guides(fill = guide_legend(reverse=TRUE)) +
  theme(
    axis.text.y = element_blank()
  )

####################### Amplicon choice
##### Amplicon choice

ampdata <- read.delim('new_analysis/evident/tech.txt', sep='\t') %>%
  select(srr, srs, project, year, worldregion, amplicon, beating) %>%
  filter(amplicon != 'None')

ampcount <- ampdata %>%
  group_by(year, amplicon) %>%
  summarise(samples=n())

# check regional numbers
regamps <- as.data.frame(table(ampdata$worldregion, ampdata$amplicon))
# fill in the blank years so they stay level for that
# year instead of disappearing
fillin <- expand.grid(
    min(ampcount$year):max(ampcount$year),
    unique(ampcount$amplicon)
  ) %>%
  mutate(dummy=1) %>%
  rename(year=Var1, amplicon=Var2) %>%
  full_join(ampcount, by=c('year','amplicon')) %>%
  select(year, amplicon, samples) %>%
  replace_na(list(samples=0))

toplot.time <- fillin %>%
  mutate(keep = amplicon %in% c('v1-v2','v3-v4','v4','v3')) %>%
  mutate(plotamplicon = ifelse(keep, amplicon, 'OTHER')) %>%
  group_by(plotamplicon, year) %>%
  summarise(
    samples=sum(samples)
  ) %>%
  arrange(year) %>%
  mutate(
    running_samples=cumsum(samples)
  )

amps <- c('v4','v3-v4','v1-v2','v3','OTHER')

ampcolors <- setNames(c(
  "#e66101",
  "#fdb863",
  "#b2abd2",
  "#5e3c99",
  "#404040"
), amps)
# colors:
ampscale <- scale_fill_manual(values=ampcolors,
                              aesthetics=c('fill', 'color'))

amps_over_time <- ggplot(toplot.time, aes(x=year, y=running_samples, color=plotamplicon)) +
  geom_line(linewidth=1) +
  labs(x='Year', y='Cumulative samples') +
  ampscale +
  scale_x_continuous(breaks=seq(2013,2021, 2), expand=c(0,0.021)) +
  scale_y_continuous(labels=scales::comma) +
  theme_bw() +
  figtheme +
  theme(
    legend.position = 'none',
    plot.margin = unit(c(1,1,1,1), "cm")
  )

total_per_amp <- ggplot(fillin, aes(x=amplicon, y=samples, fill=amplicon)) +
  geom_bar(stat='identity') +
  labs(x='Amplicon',y='Total samples') +
  ampscale +
  scale_y_continuous(labels=scales::comma) +
  theme_bw() +
  figtheme +
  theme(
    legend.position='none',
    axis.text = element_text(size=11, angle=45, hjust=1),
  )

###########
# COMPOSITION
####################

cats <- c('Bead beating*', 'World region', 'Amplicon')

alpha <- data.frame(
  cat = factor(cats, levels=rev(cats)),
  effect = c((0.43513745431541934/2), 0.23733124970037736, 0.14014893997054448)
)

plot.alpha <- ggplot(alpha, aes(y=cat, x=effect)) +
  geom_bar(stat='identity') +
  labs(x="Cohen's f", y='Factor', title='Alpha diversity') +
  scale_x_continuous(limits=c(0, 0.25)) +
  theme_bw() +
  figtheme

##### Samples per project

samplesperproject <- read.delim('new_analysis/samples_project.tsv', sep='\t')

plot.samplesproject <- ggplot(samplesperproject, aes(x=samplecount)) +
  geom_histogram(bins=30) +
  labs(x='Samples in project', y='Project count') +
  scale_x_continuous(labels=scales::comma) +
  theme_bw() +
  figtheme

median(samplesperproject$samplecount)
nrow(samplesperproject[samplesperproject$samplecount<300,])



##### Differential abundance by amplicon
ampcounts <- read.delim('new_analysis/sam/amplicon_counts.tsv') %>%
  filter(amplicon != 'None') %>%
  mutate(label = scales::comma(count)) %>%
  # add an empty row to line up with the bead beating row of the heat map
  rbind(data.frame(
    amplicon=c(''),
    count=c(0),
    label=c('')
  ))

# filter out amplicons without enough samples
keepers <- ampcounts[ampcounts$count >= 100,]$amplicon

data.amp <- read.delim('diff_abundance_results_20240705.tsv')

pretty_names <- function(names) {
  # turn periods into spaces
  names <- gsub('^[^\\.]+\\.[^\\.]+\\.[^\\.]+\\.[^\\.]+\\.[^\\.]+\\.(.+)$', '\\1', names)
  names <- gsub('NA', '(Unassigned)', names)
  names <- gsub('\\.\\.', '.', names)
  names <- gsub('\\.', ' ', names)
  names <- gsub('^ ', '', names)
}

# for each taxon, find the gap between the
# largest effect and the smallest
find_diff <- function(maxe, mine) {
  if(maxe < 0) {
    return(-mine)
  }
  if(mine > 0) {
    return(maxe)
  }
  return(maxe-mine)
}

diffs <- data.amp %>%
  filter(stat != 'uses.beating') %>%
  mutate(
    sigEstimate = if_else(adjP <= 0.05, Estimate, 0)
  ) %>%
  summarise(
    .by=taxon,
    maxe = max(Estimate),
    mine = min(Estimate),
    diff = find_diff(maxe,mine),
    maxe.sig = max(sigEstimate),
    mine.sig = min(sigEstimate),
    diff.sig = maxe.sig - mine.sig
  )
# fix values where all effects are on the same side of zero


toplot.diffs <- data.amp %>%
  mutate(
    sig = if_else(adjP <= 0.05, '*', ''),
    pretty.taxon = pretty_names(.$taxon),
    phylum = gsub('^[^\\.]+\\.([^\\.]+)\\..+$', '\\1', .$taxon),
    pretty.stat = gsub('uses\\.amp', 'v', .$stat)
  ) %>%
  mutate(
    pretty.stat = gsub('v(\\d)(\\d)', 'v\\1-v\\2', .$pretty.stat)
  ) %>%
  mutate(
    pretty.stat = gsub('uses.beating','Bead beating', .$pretty.stat)
  ) %>%
  left_join(diffs, by='taxon') %>%
  arrange(taxon) %>%
  mutate(
    pretty.taxon = factor(pretty.taxon, levels=unique(.$pretty.taxon))
    #phylum = factor(phylum, levels=unique(.$phylum))
  ) %>%
  filter(pretty.stat %in% keepers)


terp <- ggplot(toplot.diffs) +
  geom_tile(
    aes(
      x=pretty.taxon, y=pretty.stat, fill=Estimate,
      height=0.9, width=0.9
    ),
    alpha=1
  ) +
  scale_fill_gradient2(low='#2c7bb6', mid='#ffffff', high='#d7191c', midpoint=0) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  geom_text(
    aes(x=pretty.taxon, y=pretty.stat, label=sig),
    vjust=0.8
  ) +
  labs(x='Taxon', y='Amplicon') +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill='#b5b5b5'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.length.x=unit(0,'lines'),
    plot.margin = margin(1,0,-1,1,'line')
  )

plot.ampcount <- ampcounts %>%
  filter(amplicon %in% keepers) %>%  
  ggplot(aes(x=count, y=amplicon, label=label)) +
  geom_bar(stat='identity') +
  scale_x_continuous(
    breaks=c(0, 22000, max(ampcounts$count)),
    position='top',
    labels=scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  theme_bw() +
  labs(x='Samples') +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length.y=unit(0,'lines'),
    plot.margin = margin(1,1,-1,-1,'line')
  )


phylum_palette <- setNames(
  c('#BCD2EE', '#832161',     '#06D6A0',       '#E88873',   '#6153CC',       'gray', '#bc80bd',    '#fccde5',       '#1f78b4',     'gray',         '#fdb462' ),
  c('Bacillota','Pseudomonadota','Actinomycetota','Bacteroidota','Desulfobacterota','other','Euryarchaeota','Campilobacterota','Fusobacteriota','Spirochaetota','Verrucomicrobiota')
)

phylum_scale <- scale_fill_manual(values=phylum_palette, aesthetics=c('fill'))

toplot.diffs[toplot.diffs$phylum=='Firmicutes',]$phylum <- 'Bacillota'
toplot.diffs[toplot.diffs$phylum=='Proteobacteria',]$phylum <- 'Pseudomonadota'
toplot.diffs[toplot.diffs$phylum=='Actinobacteriota',]$phylum <- 'Actinomycetota'

toplot.diffs$phylum = factor(toplot.diffs$phylum, levels=unique(toplot.diffs$phylum))

bottom <- ggplot(toplot.diffs) +
  geom_tile(
    aes(x=pretty.taxon,y='Phylum',fill=phylum)
  ) +
  scale_y_discrete(expand=c(0,0)) +
  phylum_scale +
  geom_text(
    aes(x=pretty.taxon, y='Phylum', label=sig),
    vjust=0.8
  ) +
  labs(x='Taxon', fill='Phylum') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=8, angle=45, hjust=1),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(-1,-1,1,1,'line')
  )

plot.diff <- terp + plot.ampcount +
  bottom + plot_spacer() +
  plot_layout(
    heights=c(11,1),
    widths=c(11,1),
    guides = 'collect'
  ) &
  theme(
    legend.position = 'bottom'
  )

### Assemble everything
((assaycombined) | (proportional | raw)) /
  plot_spacer() /
( (amps_over_time / total_per_amp) | (plot.samplesproject / plot.alpha)) /
  plot_spacer() /
plot.diff +
plot_annotation(tag_levels = list(c('A', '', 'B','C','D','E','F','G','H'))) +
  plot_layout(heights=c(3, -0.3, 6, -0.3, 2))

ggsave('figures/fig3.pdf', device='pdf', height=13, width=13, units='in')
