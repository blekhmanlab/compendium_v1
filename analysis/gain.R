library(dplyr)
library(ggplot2)
library(patchwork)

combos <- read.delim('combinations_long.tsv')

process_iteration <- function(filename) {
  data <- read.delim(filename) %>%
    rename(entry=X.SampleID)
  
  regions <- data %>%
    slice_head(n=7)
  
  combined <- combos %>%
    inner_join(data, by=join_by(base == entry)) %>%
    rename(div.base=faith_pd) %>%
    inner_join(data, by='entry') %>%
    rename(div.combined=faith_pd) %>%
    mutate(gain = div.combined-div.base)
}
regionlist <- rev(c(
  'Europe and Northern America',
  'Australia/New Zealand',
  'Central and Southern Asia',
  'Eastern and South-Eastern Asia',
  'Oceania',
  'Northern Africa and Western Asia',
  'Latin America and the Caribbean',
  'Sub-Saharan Africa'
))


no.subsample <- process_iteration('base.tsv') %>%
  filter(base != 'Oceania') %>%
  filter(addition != 'Oceania')
no.subsample$base = factor(no.subsample$base, levels=rev(regionlist))
no.subsample$addition = factor(no.subsample$addition, levels=regionlist)

result_files <- dir(path='./gain_results', pattern="*.tsv")

recorded <- data.frame(
  base=character(),
  addition=character(),
  entry=character(),
  div.base=numeric(),
  div.combined=numeric(),
  gain=numeric(),
  run=numeric()
)

for(i in 1:length(result_files)) {
  uno <- process_iteration(paste0('gain_results/',result_files[i]))
  uno$run <- i
  recorded <- rbind(recorded, uno)
}

all_iters <- recorded %>%
  filter(base != 'Oceania') %>%
  filter(addition != 'Oceania') %>%
  mutate(base = factor(.$base, levels=rev(regionlist)),
         addition = factor(.$addition, levels=regionlist),
         prop = (100*gain) / div.base
  )

# diversity of the regions in all the iterations
plot.base <- ggplot(all_iters, aes(x=div.base, y=factor(base, levels=regionlist))) +
  geom_boxplot() +
  coord_cartesian(xlim=c(0,90000)) +
  scale_x_continuous(
    breaks=c(10000, 40000, 80000)
  ) +
  labs(
    title=paste0('Mean diversity (n=',length(result_files),')'),
    x='Faith\'s PD',x='Base',y='Addition'
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size=11)
  )

# gain!
plot.iters <- ggplot(all_iters, aes(x=base, y=addition, fill=gain)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  labs(title='Gain', y='Addition',x='', fill='Faith\'s PD') +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=11),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position='right'
  )

plot.proportional <- ggplot(all_iters, aes(x=base, y=addition, fill=prop)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  labs(title='Proportional gain', x='Base',y='Addition',fill='Gain%') +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=11),
    axis.text.x = element_text(size=11, angle=45, hjust=1),
    panel.grid = element_blank(),
    legend.position='right'
  )

plotpile <- (plot.iters / plot.proportional) /  (plot.base) +
  plot_annotation(tag_levels='A')

ggsave(plot=plotpile,
  file='gain.pdf',
  device='pdf', height=12, width=7, units='in')

# Figure out averages
all_iters.summary <- all_iters %>%
  group_by(base, addition) %>%
  summarise(
    gain=mean(gain),
    prop=mean(prop)
  )

# individual regions
region.summary <- all_iters %>%
  group_by(base) %>%
  summarise(
    faithpd=mean(div.base)
  )
