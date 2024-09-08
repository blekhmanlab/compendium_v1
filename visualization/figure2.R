library(dplyr)
library(ggplot2)
library(vegan) # for shannon
library(scales)
library(patchwork)
library(ggridges)
library(ggnewscale) # for gray backgrounds of heat maps

regionlist <- c(
  'Europe and Northern America',
  'Australia/New Zealand',
  'Central and Southern Asia',
  'Eastern and South-Eastern Asia',
  'Oceania',
  'Northern Africa and Western Asia',
  'Latin America and the Caribbean',
  'Sub-Saharan Africa'
)

regioncolors <- setNames(c(
  "#cc79a7",
  "#009e74",
  "#56b3e9",
  "#0071b2",
  "#000000",
  "#e69d00",
  "#f0e442",
  "#d55e00"
), regionlist)
# colors:
region_scale <- scale_fill_manual(values=regioncolors, aesthetics=c('fill'))

##############
# PCOA
#############
nmds <- readRDS('nmds.rds')
points <- readRDS('pcoa_points.rds')

# SCREE PLOT
scree <- data.frame(
  axis=9-row_number(nmds$eigen),
  eigen=nmds$eigen,
  eigen.rel=nmds$eigen/sum(nmds$eigen)
)
scree$cumulative <- cumsum(scree$eigen.rel)

eigenvalues <- ggplot(scree, aes(x=axis,y=eigen)) +
  geom_bar(stat='identity') +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,8), expand=c(0.01,0.01)) +
  labs(x='',
       y='Eigenvalue',
       title='Scree plot, MDS all samples'
  ) +
  theme(
    panel.grid.minor = element_blank()
  )
cumulative <- ggplot(scree, aes(x=axis,y=cumulative)) +
  geom_line(linewidth=1) +
  theme_bw() +
  scale_x_continuous(limits=c(1,8), breaks=seq(1,8), expand=c(0.065,0.065)) +
  scale_y_continuous(limits=c(0,1.05), labels = label_percent()) +
  labs(x='Component number',
       y='Variance explained (cumulative)'
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

eigen_total <- sum(nmds$eigen) / 100 # Needed for displaying % of variance explained

# Supplementary Figure 6
screeplot <- eigenvalues / cumulative
ggsave('scree.pdf', plot=screeplot, device='pdf')

######
# Annotate the samples with country/region info
######
meta <- read.delim('sample_metadata.tsv', sep='\t') %>%
  mutate(sample=paste(project, srr, sep='_'))

points.annotated <- points %>%
  left_join(
    meta,
    by='sample'
  )

points.annotated$region <- factor(points.annotated$region, levels=c(
  'Europe and Northern America',
  'Eastern and South-Eastern Asia',
  'Sub-Saharan Africa',
  'Central and Southern Asia',
  'Australia/New Zealand',
  'Northern Africa and Western Asia',
  'Latin America and the Caribbean',
  'Oceania',
  'unknown'
))

### Combined region scatter plot

# Put bigger regions toward the bottom layer so we can
# see dots from other regions
reordered.asia <- points.annotated %>%
  filter(region=='Eastern and South-Eastern Asia')

reordered.other <- points.annotated %>%
  filter(!region %in% c('unknown','Europe and Northern America','Eastern and South-Eastern Asia'))

reordered <- points.annotated %>%
  filter(region=='Europe and Northern America') %>%
  rbind(reordered.asia) %>%
  rbind(reordered.other)

xlims <- c(-10.5, 12.5)
ylims <- c(-8, 11)
regionscatter <- ggplot(reordered, aes(x=mds1, y=mds2, color=region)) +
  geom_point(size=0.1, alpha=0.5) +
  theme_bw() +
  labs(
    x='',
    y=paste('MDS2 (', round(nmds$eigen[2]/eigen_total,1), '%)', sep=''), 
    color="Region",
    title='All samples'
  ) +
  scale_x_continuous(limits=xlims, expand=c(0,0)) +
  scale_y_continuous(limits=ylims, expand=c(0,0)) +
  scale_fill_manual(values=regioncolors, aesthetics=c('color')) +
  guides(
    color = guide_legend(
      override.aes = list(size=10),
      title.position = "top"
    )
  ) +
  theme(
    legend.position='none',
    plot.title = element_text(size=10)
  )

###################
# Region-level density plots
##################
region_heat <- function(reg, xlab, ylab, newxlims, newylims, bins, v1, v2) {
  toplot <- points.annotated %>%
    filter(region == reg) %>%
    select(matches(c(v1,v2))) # only grab the two dimensions specified

  background <- points.annotated %>%
    select(matches(c(v1,v2)))

  colnames(toplot) <- c('v1','v2')
  colnames(background) <- c('v1','v2')
  heat <- ggplot(toplot, aes(x=v1, y=v2) ) +
            geom_bin_2d(data=background, bins=bins) +
            scale_colour_gradient(low='#d1d1d1', high='#d1d1d1', aesthetics='fill') +
            new_scale('fill') +
            geom_bin_2d(data=toplot, bins = 30)+
            scale_fill_continuous(type = "viridis") +
            scale_x_continuous(limits=newxlims, expand=c(0,0)) +
            scale_y_continuous(limits=newylims, expand=c(0,0)) +
            geom_vline(xintercept=0) +
            geom_hline(yintercept=0) +
            theme_bw() +
            labs(
              x=xlab, 
              y=ylab,
              title=gsub('and ', 'and\n', reg)
            ) +
            theme(
              legend.position='none',
              plot.title = element_text(size=10)
            )
  return(heat)
}

a1 <- region_heat('Europe and Northern America', '', '', xlims, ylims, 30, 'mds1','mds2')
b1 <- region_heat('Eastern and South-Eastern Asia', '', '', xlims, ylims, 30, 'mds1','mds2')
c1 <- region_heat('Sub-Saharan Africa', '', '', xlims, ylims, 30, 'mds1','mds2')
d1 <- region_heat('Central and Southern Asia', paste('MDS1 (', round(nmds$eigen[1]/eigen_total,1), '%)', sep=''),'MDS2', xlims, ylims, 30, 'mds1','mds2')
e1 <- region_heat('Australia/New Zealand', '', '', xlims, ylims, 30, 'mds1','mds2')
f1 <- region_heat('Northern Africa and Western Asia', '', '', xlims, ylims, 30, 'mds1','mds2')
g1 <- region_heat('Latin America and the Caribbean', '', '', xlims, ylims, 30, 'mds1','mds2')

########### Supplemental   ############
build_page <- function(v1, v2) {
  newbincount <- 50
  # we need the index of the axis being described so 
  # we can get the variance explained
  i1 <- strtoi(substr(v1, 4, 4))
  i2 <- strtoi(substr(v2, 4, 4))

  findmax <- points.annotated %>%
    select(matches(c(v1,v2))) # only grab the two dimensions specified
  colnames(findmax) <- c('v1','v2')
  subxlims <- c(min(findmax$v1)*1.1, max(findmax$v1)*1.1)
  subylims <- c(min(findmax$v2)*1.1, max(findmax$v2)*1.1)
  print(subxlims)
  print(subylims)
  a1 <- region_heat('Europe and Northern America',
    paste(toupper(v1), ' (', round(nmds$eigen[i1]/eigen_total,1), '%)', sep=''),
    paste(toupper(v2), ' (', round(nmds$eigen[i2]/eigen_total,1), '%)', sep=''),
    subxlims,subylims, newbincount,
    v1, v2
  )
  b1 <- region_heat('Eastern and South-Eastern Asia', '', '',subxlims,subylims,newbincount,v1, v2)
  c1 <- region_heat('Sub-Saharan Africa', '', '',subxlims,subylims,newbincount,v1, v2)
  d1 <- region_heat('Central and Southern Asia', '', '',subxlims,subylims,newbincount, v1, v2)
  e1 <- region_heat('Australia/New Zealand', '', '',subxlims,subylims,newbincount,v1, v2)
  f1 <- region_heat('Northern Africa and Western Asia', '', '',subxlims,subylims,newbincount,v1, v2)
  g1 <- region_heat('Latin America and the Caribbean', '', '',subxlims,subylims,newbincount,v1, v2)

  return(a1+b1+c1+d1+e1+f1+g1)
}

page1 <- build_page('mds1','mds2')
page2 <- build_page('mds1','mds3')
page3 <- build_page('mds2','mds3')
page4 <- build_page('mds1','mds4')
page5 <- build_page('mds2','mds4')
page6 <- build_page('mds3','mds4')

ggsave(plot=page1, file='S8_pcoa.pdf', device='pdf', height=9, width=8, units='in')
ggsave(plot=page2, file='S9_pcoa.pdf', device='pdf', height=9, width=8, units='in')
ggsave(plot=page3, file='S10_pcoa.pdf', device='pdf', height=9, width=8, units='in')
ggsave(plot=page4, file='S11_pcoa.pdf', device='pdf', height=9, width=8, units='in')
ggsave(plot=page5, file='S12_pcoa.pdf', device='pdf', height=9, width=8, units='in')
ggsave(plot=page6, file='S13_pcoa.pdf', device='pdf', height=9, width=8, units='in')

##############################


mds_ridges <- function(mds, ylab) {
  toplot <- points.annotated %>%
    select(sample, region, mds) %>%
    rename(all_of(c(selected=mds))) %>%
    filter(!region %in% c('', 'unknown', 'Oceania'))
  
  toplot$region <- factor(toplot$region,
                          levels=rev(c(
                            'Europe and Northern America',
                            'Eastern and South-Eastern Asia',
                            'Sub-Saharan Africa',
                            'Central and Southern Asia',
                            'Australia/New Zealand',
                            'Northern Africa and Western Asia',
                            'Latin America and the Caribbean',
                            'Oceania'
                          )))

  ridgeplot <- ggplot(toplot, aes(x=selected, y=region, fill=region)) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
    theme_bw() +
    region_scale +
    labs(x=toupper(mds),y='')
  
  if(ylab) {
    return(ridgeplot +
             theme(
               legend.position='none',
               plot.title = element_text(size=10),
               axis.text = element_text(size=11),
               rect = element_rect(fill='transparent')
             )
    )
  }
  return(ridgeplot +
           theme(
             legend.position='none',
             plot.title = element_text(size=10),
             axis.text.x = element_text(size=11),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             rect = element_rect(fill='transparent')
           )
  )
}

# Figure 2F
pmds1 <- mds_ridges('mds1', T)
pmds2 <- mds_ridges('mds2', FALSE)
pmds3 <- mds_ridges('mds3', FALSE)
pmds4 <- mds_ridges('mds4', FALSE)
# (for supplement)
pmds5 <- mds_ridges('mds5', T)
pmds6 <- mds_ridges('mds6', FALSE)
pmds7 <- mds_ridges('mds7', FALSE)
pmds8 <- mds_ridges('mds8', FALSE)

allmds <- pmds1 + pmds2 + pmds3 + pmds4 + 
  pmds5 + pmds6 + pmds7 + pmds8 +
  plot_layout(ncol=4, nrow=2)
ggsave('allmds.pdf', plot=allmds, device='pdf', height=5.6, width=10, units='in')


# Test differences in mds distributions
if(FALSE) {
  kstest_region <- function(region, mds) {
    ks.test(
      points.annotated[points.annotated$region==region,][[mds]],
      points.annotated[!points.annotated$region==region,][[mds]]
    )
  }
  kstest_region('Europe and Northern America', 'mds1')
  kstest_region('Europe and Northern America', 'mds2')
  kstest_region('Europe and Northern America', 'mds3')
  kstest_region('Europe and Northern America', 'mds4')
  
  kstest_region('Latin America and the Caribbean', 'mds1')
  kstest_region('Latin America and the Caribbean', 'mds2')
  kstest_region('Latin America and the Caribbean', 'mds3')
  kstest_region('Latin America and the Caribbean', 'mds4')
  
  kstest_region('Australia/New Zealand', 'mds1')
  
  kstest_regions <- function(region1, region2, mds) {
    ks.test(
      points.annotated[points.annotated$region==region1,][[mds]],
      points.annotated[points.annotated$region==region2,][[mds]]
    )
  }
  kstest_regions('Europe and Northern America', 'Australia/New Zealand', 'mds1')
  kstest_regions('Sub-Saharan Africa', 'Eastern and South-Eastern Asia', 'mds4')
}

########## Read depth and diversity
working <- readRDS('filtered.rds')

meta$join <- paste(meta$project, meta$srr, sep='_')

totals <- as.data.frame(rowSums(working))
totals$sample <- rownames(totals)
totals <- rename(totals, total='rowSums(working)')

data <- meta %>%
  select(join, region) %>%
  inner_join(totals, by=join_by(join==sample))

data$region <- factor(data$region,
                      levels=rev(c(
                        'Europe and Northern America',
                        'Eastern and South-Eastern Asia',
                        'Sub-Saharan Africa',
                        'Central and Southern Asia',
                        'Australia/New Zealand',
                        'Northern Africa and Western Asia',
                        'Latin America and the Caribbean',
                        'Oceania',
                        'unknown'
                      )))

# Figure 2D
depth <- data  %>%
  filter(!region %in% c('unknown','','Oceania')) %>%
  ggplot(aes(x=total, y=region, fill=region)) +
  geom_violin(draw_quantiles=0.5) +
  scale_x_log10(breaks=c(1e4, 1e5, 1e6, 1e7), labels=scales::label_log()) +
  labs(x='Read depth',y='') +
  region_scale +
  theme_bw() +
  theme(
    legend.position='none',
    strip.background = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    rect = element_rect(
      fill='transparent'
    ),
    plot.tag.position = c(0.1,0.965)
  )

####################
# (Bar plot: samples per region)
###################
rename_region <- function(data, oldname, newname) {
  data[data$printname==oldname,]$printname = newname
  return(data)
}

samplecounts <- data %>%
  filter(!region %in% c('Oceania','unknown'))

samplecounts <- table(samplecounts$region) %>%
  as.data.frame() %>%
  filter(!Var1 %in% c('Oceania','unknown')) %>%
  mutate(printname=as.character(Var1))

samplecounts$Var1 <- factor(samplecounts$Var1, levels=rev(c(
  'Europe and Northern America',
  'Eastern and South-Eastern Asia',
  'Sub-Saharan Africa',
  'Central and Southern Asia',
  'Australia/New Zealand',
  'Northern Africa and Western Asia',
  'Latin America and the Caribbean',
  'Oceania'
)))

# Shorten names for display
samplecounts <- samplecounts %>%
  rename_region('Europe and Northern America', 'Europe and N. America') %>%
  rename_region('Eastern and South-Eastern Asia', 'Eastern and SE Asia') %>%
  rename_region('Central and Southern Asia', 'Central and S. Asia') %>%
  rename_region('Northern Africa and Western Asia', 'N. Africa and W. Asia') %>%
  rename_region('Latin America and the Caribbean', 'Latin America/Caribbean')

samplecount_plot <- ggplot(samplecounts, aes(y=Var1, x=Freq, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.5) +
  geom_text(
    aes(
      x=1000,
      y=Var1,
      #y=c(rep(1000,9), rep(30000, 1)),
      label=printname
    ),
    hjust=0, size=4, color='black'
  ) +
  scale_x_continuous(
    breaks=c(0, 2e4, 5e4, 1e5),
    limits = c(0,101000),
    expand=c(0,5000),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  region_scale +
  theme_bw() +
  labs(x='Samples',y='Region') +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=11),
    rect = element_rect(fill='transparent')
  )

##########
# Diversity
#########
# make a version of the dataset with only samples that
# have a regional annotation
working$sample <- rownames(working)
only_regional <- working %>%
  inner_join(meta, by=join_by(sample==join)) %>%
  filter(!region %in% c('unknown','','Oceania')) %>%
  select(!any_of(colnames(meta)))
rownames(only_regional) <- only_regional$sample
only_regional$sample <- NULL


only_regional.annotated <- only_regional %>%
  mutate(join=rownames(only_regional)) %>%
  inner_join(data, by='join')
rownames(only_regional.annotated) <- only_regional.annotated$join

unrarefied <- only_regional.annotated %>%
  select(!c('join','region','total','sample.y')) %>%
  diversity(index='shannon') %>%
  as.data.frame() %>%
  rename(shannon='.') %>%
  mutate(region=only_regional.annotated$region)

unrarefied$region <- factor(unrarefied$region,
                            levels=rev(c(
                              'Europe and Northern America',
                              'Eastern and South-Eastern Asia',
                              'Sub-Saharan Africa',
                              'Central and Southern Asia',
                              'Australia/New Zealand',
                              'Northern Africa and Western Asia',
                              'Latin America and the Caribbean',
                              'Oceania'
                            )))

###### Supplemental figure for rarefaction results
# rarefaction_diversity.rds is calculated by diversity_rarefaction.R

# We do this here because panel 2C requires the "alphadata" variable
alpha <- readRDS('rarefaction_diversity.rds') %>%
  filter(!region %in% c('unknown','','Oceania'))

alpha$region <- factor(alpha$region,
                        levels=rev(c(
                          'Europe and Northern America',
                          'Eastern and South-Eastern Asia',
                          'Sub-Saharan Africa',
                          'Central and Southern Asia',
                          'Australia/New Zealand',
                          'Northern Africa and Western Asia',
                          'Latin America and the Caribbean',
                          'Oceania'
                        )))

conf_interval <- function(data, side) {
  result <- t.test(data)
  
  if(side == 'bottom') {
    # bottom of the interval
    return(result$conf.int[1])
  }
  # otherwise, return the top
  return(result$conf.int[2])
}

alphadata <- alpha %>%
  group_by(region) %>%
  summarise(
    sd=sd(diversity),
    mean=mean(diversity),
    n=length(diversity),
    median=median(diversity),
    conf_bottom = conf_interval(diversity, 'bottom'),
    conf_top = conf_interval(diversity, 'top'),
  ) %>%
  mutate(
    se = sd / sqrt(n)
  )

# Supplementary table
write.table(alphadata, file='alpha_summary.tsv',
            sep='\t', row.names = F)

rareplot <- ggplot(alphadata,
         aes(x=mean, y=region, xmin=conf_bottom,
             xmax=conf_top, color=region)) +
  geom_point() +
  geom_errorbar() +
  scale_fill_manual(values=regioncolors, aesthetics=c('color')) +
  theme_bw() +
  labs(x='Shannon diversity', y='Region') +
  theme(
    legend.position='none'
  )
ggsave('rarefaction.pdf', plot=rareplot, device='pdf', height=10.2, width=9.6, units='in')



depthmedians <- data %>%
  group_by(region) %>%
  summarise(median=median(total))
########################################

rawdiv <- ggplot(unrarefied, aes(x=shannon, y=region, fill=region)) +
  geom_violin(draw_quantiles = 0.5) +
  region_scale +
  geom_point(data=alphadata, mapping=aes(y=region, x=mean), size=2) +
  theme_bw() +
  theme(
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    rect = element_rect(fill='transparent')
  ) +
  labs(x='Shannon diversity')

# Supplemental figure for other alpha diversity measures
unrarefied.simpson <- only_regional.annotated %>%
  select(!c('join','region','total','sample.y')) %>%
  diversity(index='simpson') %>%
  as.data.frame() %>%
  rename(simpson='.') %>%
  mutate(region=only_regional.annotated$region)

unrarefied.simpson$region <- factor(unrarefied.simpson$region,
                                    levels=rev(c(
                                      'Europe and Northern America',
                                      'Eastern and South-Eastern Asia',
                                      'Sub-Saharan Africa',
                                      'Central and Southern Asia',
                                      'Australia/New Zealand',
                                      'Northern Africa and Western Asia',
                                      'Latin America and the Caribbean',
                                      'Oceania'
                                    )))
#####
unrarefied.specnumber <- only_regional.annotated %>%
  select(!c('join','region','total','sample.y')) %>%
  specnumber() %>%
  as.data.frame() %>%
  rename(specnumber='.') %>%
  mutate(region=only_regional.annotated$region)

unrarefied.specnumber$region <- factor(unrarefied.specnumber$region,
                                       levels=rev(c(
                                         'Europe and Northern America',
                                         'Eastern and South-Eastern Asia',
                                         'Sub-Saharan Africa',
                                         'Central and Southern Asia',
                                         'Australia/New Zealand',
                                         'Northern Africa and Western Asia',
                                         'Latin America and the Caribbean',
                                         'Oceania'
                                       )))
#####
rawdiv.shannon <- ggplot(unrarefied, aes(x=shannon, y=region, fill=region)) +
  geom_boxplot(outlier.shape = NA) +
  region_scale +
  theme_bw() +
  theme(
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=11),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    rect = element_rect(fill='transparent')
  ) +
  labs(x='Shannon diversity')

rawdiv.simpson <- ggplot(unrarefied.simpson, aes(x=simpson, y=region, fill=region)) +
  geom_boxplot(outlier.shape = NA) +
  region_scale +
  theme_bw() +
  theme(
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=11),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  ) +
  labs(x='Simpson diversity')

rawdiv.specnumber <- ggplot(unrarefied.specnumber, aes(x=specnumber, y=region, fill=region)) +
  geom_boxplot(outlier.shape = NA) +
  region_scale +
  coord_cartesian(xlim=c(0,180)) +
  theme_bw() +
  theme(
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=11),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  ) +
  labs(x='Species count')

divplot <- rawdiv.shannon / rawdiv.simpson / rawdiv.specnumber

ggsave('diversity.pdf', plot=divplot, device='pdf', height=10.2, width=9.6, units='in')




###############
# Figure 2A
###############

library(rnaturalearth)
library(sf)
library(grid)

regions <- read.csv('regions.csv', header=F) %>%
  rename(iso_a2=V1, country=V2, region=V3)
regions[regions$country=='Namibia',]$iso_a2 = 'NA'

toplot <- ne_countries(scale='medium', type='map_units', returnclass = 'sf') %>%
  left_join(regions, by = c('iso_a2'))

# How many entries don't have a regional assignment?
filter(toplot, is.na(region)) %>% nrow()
# Clean up subregions that don't line up exactly with our data
table(toplot[is.na(toplot$region),]$subregion)

toplot[is.na(toplot$region) & toplot$subregion=='Caribbean',]$region <- 'Latin America and the Caribbean'
toplot[is.na(toplot$region) & toplot$subregion=='Eastern Africa',]$region <- 'Sub-Saharan Africa'
toplot[is.na(toplot$region) & toplot$subregion=='Northern Europe',]$region <- 'Europe and Northern America'
toplot[is.na(toplot$region) & toplot$subregion=='Polynesia',]$region <- 'Oceania'
toplot[is.na(toplot$region) & toplot$subregion=='South-Eastern Asia',]$region <- 'Eastern and South-Eastern Asia'
toplot[is.na(toplot$region) & toplot$subregion=='Western Asia',]$region <- 'Northern Africa and Western Asia'
toplot[is.na(toplot$region) & toplot$subregion=='Australia and New Zealand',]$region <- 'Australia/New Zealand'
toplot[is.na(toplot$region) & toplot$subregion=='Melanesia',]$region <- 'Oceania'
toplot[is.na(toplot$region) & toplot$subregion=='Western Europe',]$region <- 'Europe and Northern America'
toplot[is.na(toplot$region) & toplot$subregion=='Southern Europe',]$region <- 'Europe and Northern America'
toplot[is.na(toplot$region) & toplot$subregion=='Southern Asia',]$region <- 'Central and Southern Asia'

# All fixed?
table(toplot[is.na(toplot$region),]$subregion)

toplot$region <- factor(toplot$region,
                        levels=c(
                          'Europe and Northern America',
                          'Australia/New Zealand',
                          'Central and Southern Asia',
                          'Eastern and South-Eastern Asia',
                          'Oceania',
                          'Northern Africa and Western Asia',
                          'Latin America and the Caribbean',
                          'Sub-Saharan Africa'
                        ))
mapplot <- ggplot(data=toplot) +
  geom_sf(aes(fill=region), color=NA, size=0.01) +
  coord_sf(crs = "+proj=eqearth +wktext") + # changes the projection
  labs(fill='Region membership') +
  region_scale +
  theme(
    legend.position='none',
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  )


###############
######## Final assembly
###############

toprow <- plot_spacer() + 
  samplecount_plot + plot_spacer() +
  rawdiv + plot_spacer() +
  depth + plot_spacer() +
  plot_layout(
    nrow=1,
    widths=c(9.5, 7,-0.7, 4,-0.85, 4, 0.0)
  )

# The regionheats layout looks so strange because
# we want to minimize the whitespace between plots
hgap <- -0.1
vgap <- -0.29

regionheats <- regionscatter + plot_spacer() + a1 + plot_spacer() + b1 + plot_spacer() + c1 +
  plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + 
  d1 + plot_spacer() + e1 + plot_spacer() + f1 + plot_spacer() + g1 +
  plot_layout(
    nrow=3,
    widths=c(1, hgap, 1, hgap, 1, hgap, 1),
    heights=c(1, vgap, 1)
  )

hgap2 <- -0.25
bottomrow <- pmds1 + plot_spacer() +
  pmds2 + plot_spacer() +
  pmds3 + plot_spacer() +
  pmds4 +
  plot_layout(
    nrow=1,
    widths=c(1,hgap2, 1,hgap2, 1,hgap2, 1)
  )

full <- toprow /
  wrap_elements(plot=regionheats) /
  plot_spacer() /
  wrap_elements(bottomrow) +
  plot_layout(heights=c(3,5, -0.55, 4))

ggsave('fig2.pdf', plot=full, device='pdf', height=10.2, width=8.7, units='in')
