library(ggplot2)
library(dplyr)
library(tidyr) # for pivot_longer
library(patchwork)
library(scales) # for percentages

###### PCA plots ################################

batch1 <- c('usa','china','denmark','japan','canada','united kingdom')
batch2 <- c('finland','germany','bangladesh','india','new zealand','australia')
batch3 <- c('nigeria','israel','malawi','mexico','colombia','turkey')

colors <- c(
  rgb(207, 164, 131, maxColorValue=255),
  rgb(243, 127, 61, maxColorValue=255),
  rgb(13, 206, 235, maxColorValue=255),
  rgb(219, 186, 247, maxColorValue=255),
  rgb(67, 100, 206, maxColorValue=255),
  rgb(223, 224, 59, maxColorValue=255),
  rgb(241, 51, 215, maxColorValue=255),
  rgb(179, 233, 89, maxColorValue=255),
  rgb(8, 177, 84, maxColorValue=255),
  rgb(249, 186, 206, maxColorValue=255),
  rgb(252, 220, 65, maxColorValue=255),
  rgb(237, 19, 76, maxColorValue=255),
  rgb(133, 12, 19, maxColorValue=255),
  rgb(155, 208, 204, maxColorValue=255),
  rgb(217, 193, 171, maxColorValue=255),
  rgb(250, 245, 200, maxColorValue=255),
  rgb(152, 28, 171, maxColorValue=255),
  rgb(155, 247, 194, maxColorValue=255)
)
countrycolors <- setNames(colors, c(batch1, batch2, batch3))
country_scale <- scale_fill_manual(values=countrycolors, aesthetics=c('color','fill'))


meta <- read.csv('compendium_metadata.csv') %>%
  dplyr::select(!c('X'))

cdata <- read.csv('compendium_pca.csv') %>%
  dplyr::select(c('PC1'='X0', 'PC2'='X1','PC3'='X2','PC4'='X3','PC5'='X4','PC6'='X5'))

cdata$country <- meta$country

library(ggnewscale)

# original scatter plots, working
countryplot <- function(cname) {

  batch1 <- c('usa','china','denmark','japan','canada','united kingdom')
  batch2 <- c('finland','germany','bangladesh','india','new zealand','australia')
  batch3 <- c('nigeria','israel','malawi','mexico','colombia','turkey')
  
  toplot <- cdata %>%
    filter(country == cname)
  
  col1 <- ggplot(toplot, aes(x=PC1,y=PC2, color=country)) +
    geom_point(size=0.4) +
    country_scale +
    scale_x_continuous(limits = c(min(cdata$PC1), max(cdata$PC1))) +
    scale_y_continuous(limits = c(min(cdata$PC2), max(cdata$PC2))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  col2 <- ggplot(toplot, aes(x=PC3,y=PC4, color=country)) +
    geom_point(size=0.4) +
    country_scale +
    scale_x_continuous(limits = c(min(cdata$PC3), max(cdata$PC3))) +
    scale_y_continuous(limits = c(min(cdata$PC4), max(cdata$PC4))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  col3 <- ggplot(toplot, aes(x=PC5,y=PC6, color=country)) +
    geom_point(size=0.4) +
    country_scale +
    scale_x_continuous(limits = c(min(cdata$PC5), max(cdata$PC5))) +
    scale_y_continuous(limits = c(min(cdata$PC6), max(cdata$PC6))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  
  assembled <- (col1 | col2 | col3) +
    plot_annotation(
      title = cname,
      theme = theme(
        plot.title = element_text(size = 16, hjust=0.5, face='bold')
      )
    )
  return(wrap_elements(assembled))
}


AXIS_EXPAND = 1.05 # we should make the axes slightly longer
                   # than necessary to leave room for the hexagons
# scatter with background
countryplot <- function(cname) {
  prettyname <- c(
    'usa' = 'United States',
    'china' = 'China',
    'denmark' = 'Denmark',
    'japan' = 'Finland',
    'canada' = 'Canada',
    'united kingdom' = 'United Kingdom',
    'finland' = 'Finland',
    'germany' = 'Germany',
    'bangladesh' = 'Bangladesh',
    'india' = 'India',
    'new zealand' = 'New Zealand',
    'australia' = 'Australia',
    'nigeria' = 'Nigeria',
    'israel' = 'Israel',
    'malawi' = 'Malawi',
    'mexico' = 'Mexico',
    'colombia' = 'Colombia',
    'turkey' = 'Turkey'
  )

  toplot <- cdata %>%
    filter(country == cname)
  
  col1 <- ggplot(cdata, aes(x=PC1,y=PC2)) +
    geom_hex(bins=30) +
    scale_colour_gradient(low='#878686', high='#878686', aesthetics='fill') +
    new_scale('fill') +
    geom_point(data=toplot, mapping=aes(x=PC1,y=PC2, color=country), size=0.8)+
    country_scale +
    scale_x_continuous(limits = c((min(cdata$PC1)*AXIS_EXPAND), (max(cdata$PC1)*AXIS_EXPAND))) +
    scale_y_continuous(limits = c((min(cdata$PC2)*AXIS_EXPAND), (max(cdata$PC2)*AXIS_EXPAND))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  col2 <- ggplot(cdata, aes(x=PC3,y=PC4)) +
    geom_hex(bins=30) +
    scale_colour_gradient(low='#878686', high='#878686', aesthetics='fill') +
    new_scale('fill') +
    geom_point(data=toplot, mapping=aes(x=PC3,y=PC4, color=country), size=0.8)+
    country_scale +
    scale_x_continuous(limits = c((min(cdata$PC3)*AXIS_EXPAND), (max(cdata$PC3)*AXIS_EXPAND))) +
    scale_y_continuous(limits = c((min(cdata$PC4)*AXIS_EXPAND), (max(cdata$PC4)*AXIS_EXPAND))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  col3 <- ggplot(cdata, aes(x=PC5,y=PC6)) +
    geom_hex(bins=30) +
    scale_colour_gradient(low='#878686', high='#878686', aesthetics='fill') +
    new_scale('fill') +
    geom_point(data=toplot, mapping=aes(x=PC5,y=PC6, color=country), size=0.8)+
    country_scale +
    scale_x_continuous(limits = c((min(cdata$PC5)*AXIS_EXPAND), (max(cdata$PC5)*AXIS_EXPAND))) +
    scale_y_continuous(limits = c((min(cdata$PC6)*AXIS_EXPAND), (max(cdata$PC6)*AXIS_EXPAND))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  
  assembled <- (plot_spacer() | col1 | col2 | col3 | plot_spacer()) +
    plot_layout(nrow=1, widths=c(-0.2, 1,1,1, -0.2)) +
    plot_annotation(
      title = prettyname[cname],
      theme = theme(
        plot.title = element_text(size = 16, hjust=0.5)#, face='bold')
      )
    )
  return(wrap_elements(assembled))
}

# heatmaps
countryplot <- function(cname) {
  prettyname <- c(
    'usa' = 'United States',
    'china' = 'China',
    'denmark' = 'Denmark',
    'japan' = 'Finland',
    'canada' = 'Canada',
    'united kingdom' = 'United Kingdom',
    'finland' = 'Finland',
    'germany' = 'Germany',
    'bangladesh' = 'Bangladesh',
    'india' = 'India',
    'new zealand' = 'New Zealand',
    'australia' = 'Australia',
    'nigeria' = 'Nigeria',
    'israel' = 'Israel',
    'malawi' = 'Malawi',
    'mexico' = 'Mexico',
    'colombia' = 'Colombia',
    'turkey' = 'Turkey'
  )

  toplot <- cdata %>%
    filter(country == cname)
  
  col1 <- ggplot(cdata, aes(x=PC1,y=PC2)) +
    geom_hex(bins=30) +
    scale_colour_gradient(low='#d1d1d1', high='#d1d1d1', aesthetics='fill') +
    new_scale('fill') +
    geom_hex(data=toplot, bins = 30)+
    scale_fill_continuous(type = "viridis") +
    scale_x_continuous(limits = c((min(cdata$PC1)*AXIS_EXPAND), (max(cdata$PC1)*AXIS_EXPAND))) +
    scale_y_continuous(limits = c((min(cdata$PC2)*AXIS_EXPAND), (max(cdata$PC2)*AXIS_EXPAND))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  col2 <- ggplot(cdata, aes(x=PC3,y=PC4)) +
    geom_hex(bins=30) +
    scale_colour_gradient(low='#d1d1d1', high='#d1d1d1', aesthetics='fill') +
    new_scale('fill') +
    geom_hex(data=toplot, bins = 30)+
    scale_fill_continuous(type = "viridis") +
    scale_x_continuous(limits = c((min(cdata$PC3)*AXIS_EXPAND), (max(cdata$PC3)*AXIS_EXPAND))) +
    scale_y_continuous(limits = c((min(cdata$PC4)*AXIS_EXPAND), (max(cdata$PC4)*AXIS_EXPAND))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  col3 <- ggplot(cdata, aes(x=PC5,y=PC6)) +
    geom_hex(bins=30) +
    scale_colour_gradient(low='#d1d1d1', high='#d1d1d1', aesthetics='fill') +
    new_scale('fill') +
    geom_hex(data=toplot, bins = 30)+
    scale_fill_continuous(type = "viridis") +
    scale_x_continuous(limits = c((min(cdata$PC5)*AXIS_EXPAND), (max(cdata$PC5)*AXIS_EXPAND))) +
    scale_y_continuous(limits = c((min(cdata$PC6)*AXIS_EXPAND), (max(cdata$PC6)*AXIS_EXPAND))) +
    theme_bw() +
    theme(
      legend.position='none',
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  
  assembled <- (plot_spacer() | col1 | col2 | col3 | plot_spacer()) +
    plot_layout(nrow=1, widths=c(-0.2, 1,1,1, -0.2)) +
    plot_annotation(
      title = prettyname[cname],
      theme = theme(
        plot.title = element_text(size = 16, hjust=0.5)
      )
    )
  return(wrap_elements(assembled))
}

batchplot <- function(countries) {
  hplot <- 1   # relative height of plots
  hspace <- -0.2 # gap between rows of the plot
  assembled <- countryplot(countries[1]) /
    plot_spacer() /
    countryplot(countries[2]) /
    plot_spacer() /
    countryplot(countries[3]) /
    plot_spacer() /
    countryplot(countries[4]) /
    plot_spacer() /
    countryplot(countries[5]) /
    plot_spacer() /
    countryplot(countries[6]) +
    plot_layout(heights=c(rep(c(hplot, hspace),5), hplot))

  return(assembled)
}

batch1.plot <- batchplot(batch1)
batch2.plot <- batchplot(batch2)
batch3.plot <- batchplot(batch3)

tophalf <- (batch1.plot | batch2.plot | batch3.plot)

###### Bottom row ################################
regions <- c(
  'Australia/New Zealand',
  'Central and Southern Asia',
  'Northern Africa and Western Asia',
  'Sub-Saharan Africa',
  'Latin America and the Caribbean',
  'Eastern and South-Eastern Asia',
  'Europe and Northern America'
)

regioncolors <- setNames(c(
  "#009e74",
  "#56b3e9",
  "#e69d00",
  "#d55e00",
  "#f0e442",
  "#0071b2",
  "#cc79a7"
), regions)

# colors:
region_scale <- scale_fill_manual(values=regioncolors, aesthetics=c('fill','color'))

auc <- data.frame(
  region = factor(regions, levels=rev(regions)),
  area = c(0.944, 0.942, 0.892, 0.875, 0.852, 0.831, 0.797)
)

curvedata <- read.csv('new_analysis/ashwin/auc_curves_ova_classifiers.csv') %>%
  mutate(region=factor(region, levels=rev(regions)))

bigcurve <- ggplot(curvedata, aes(x=fpr, y=tpr, color=region)) +
  geom_point() +
  geom_line() +
  region_scale +
  labs(title='Regional classifier', x='False positive rate', y='True positive rate') +
  theme_bw() +
  theme(
    legend.position='none'
  )

aucplot <- ggplot(auc, aes(x=area, y=region, fill=region)) +
  geom_bar(stat='identity', alpha=0.9) +
  region_scale +
  scale_x_continuous(expand=c(0, 0, 0.1, 0)) + # this format makes no sense?
  geom_text(aes(x=0.02, y=region, label=region), hjust=0) +
  labs(x='AUC', y='') +
  theme_bw() +
  theme(
    legend.position='none',
    axis.text.y = element_blank()
  )


results <- readRDS('country_cluster_bootstrap.100min.rds')
real <- readRDS('country_cluster_bootstrap.100min.REAL.rds')

iters <- ggplot(results, aes(x=results)) +
  geom_histogram(bins=75) +
  geom_vline(xintercept=real$DB, color='red', linewidth=1) +
  annotate(geom='text', x=real$DB+50, y=600, label=paste0('True value: ',round(real$DB, digits=3))) +
  annotate(geom='text', x=170, y=900, label=paste0('Median: ',round(median(results$results), digits=1))) +
  labs(x='Davies-Bouldin Index', y='Iterations', title='Country-level clustering') +
  theme_bw()

# "iters" is from the country_cluster_evaluation.R

bottomhalf <- (iters | bigcurve | aucplot)

# no idea why free(wrap_elements()) fixes the alignment issue
assembled <- free(wrap_elements(tophalf)) + bottomhalf +
  plot_layout(ncol=1, nrow=2, heights=c(5,1)) +
  plot_annotation(tag_levels='A') &
  theme(
    plot.tag=element_text(face='bold')
  )

ggsave(assembled, file='figures/fig6.pdf', device='pdf',
  height=15, width=12, units='in'
)