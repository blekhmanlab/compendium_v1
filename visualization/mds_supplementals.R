BACKGROUND <- reordered %>%
  filter(!region %in% c('unknown', 'Oceania'))

bincount <- 50

MAKE_PLOT <- function(mdsx, mdsy) {
  first <- paste0('mds',mdsx)
  second <- paste0('mds',mdsy)
  
  regionscatter <- ggplot(BACKGROUND, aes(
      x=eval(as.symbol(first)),
      y=eval(as.symbol(second)),
      color=region
    )) +
    geom_point(size=0.1, alpha=0.5) +
    theme_bw() +
    labs(
      x=paste(toupper(first), ' (', round(nmds$eigen[mdsx]/eigen_total,1), '%)', sep=''), 
      y=paste(toupper(second), ' (', round(nmds$eigen[mdsy]/eigen_total,1), '%)', sep=''), 
      color="Region",
      title='All samples'
    ) +
    #scale_x_continuous(limits=xlims, expand=c(0,0)) +
    #scale_y_continuous(limits=ylims, expand=c(0,0)) +
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
  
  region_heat <- function(reg) {
    toplot <- BACKGROUND %>% filter(region == reg)
    
    heat <- ggplot(toplot, aes(
        x=eval(as.symbol(first)),
        y=eval(as.symbol(second))
      )) +
      geom_bin_2d(data=BACKGROUND, bins=bincount) +
      scale_colour_gradient(low='#d1d1d1', high='#d1d1d1', aesthetics='fill') +
      new_scale('fill') +
      geom_bin_2d(data=toplot, bins = bincount)+
      scale_fill_continuous(type = "viridis") +
      #scale_x_continuous(limits=xlims, expand=c(0,0)) +
      #scale_y_continuous(limits=ylims, expand=c(0,0)) +
      geom_vline(xintercept=0) +
      geom_hline(yintercept=0) +
      theme_bw() +
      labs(
        x='', 
        y='',
        title=gsub('and ', 'and\n', reg)
      ) +
      theme(
        legend.position='none',
        plot.title = element_text(size=10)
      )
    return(heat)
  }
  a1 <- region_heat('Europe and Northern America')
  b1 <- region_heat('Eastern and South-Eastern Asia')
  c1 <- region_heat('Sub-Saharan Africa')
  d1 <- region_heat('Central and Southern Asia')
  e1 <- region_heat('Australia/New Zealand')
  f1 <- region_heat('Northern Africa and Western Asia')
  g1 <- region_heat('Latin America and the Caribbean')
  
  built <- regionscatter + a1 + b1 + c1 + d1 + e1 + f1 + g1 +
    plot_layout(nrow=2) +
    plot_annotation(paste0(toupper(first), ' vs ', toupper(second)))
  
  return(built)
}

plot12 <- MAKE_PLOT(1, 2)
plot13 <- MAKE_PLOT(1, 3)
plot23 <- MAKE_PLOT(2, 3)
plot14 <- MAKE_PLOT(1, 4)
plot24 <- MAKE_PLOT(2, 4)
plot34 <- MAKE_PLOT(3, 4)

plot12 / plot13 / plot23 / plot14 / plot24 / plot34
