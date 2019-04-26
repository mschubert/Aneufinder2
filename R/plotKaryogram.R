plotKaryogram <- function(model) {
  bins                  <- model$bins                                                                                  # GRanges with bin counts
  bins.split            <- split(bins, seqnames(bins))                                                                 # GRanges list (every item is one chromosome)
  bins.split            <- bins.split[lengths(bins.split) > 0]                                                         #   Ref can have more chr names than 1-22, X and Y so remove others
  custom.xlim           <- stats::quantile(bins$counts, 0.999)                                                         # Maximum count (xlim for every chromosome)
  if(length(custom.xlim) == 0 | is.na(custom.xlim) | is.infinite(custom.xlim)){                                        #   New code replacing function: 'get_rightxlim'
    custom.xlim         <- 1
  }
  empty_theme           <- theme(axis.line=element_blank(),axis.text=element_blank(),
                                 axis.ticks=element_blank(),axis.title.x=element_text(size=13),                        # 'axis.title.x': Text size chromosomes x-axis
                                 axis.title.y=element_blank(),legend.position="none",
                                 panel.background=element_blank(),panel.border=element_blank(),                        # Are these panel settings necessary? I don't see a change when removing them
                                 panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                                 plot.background=element_blank())
  ggplts                <- list()                                                                                      # Make plot object for every chromosome (every chromosome is one plot)
  for(chrom in names(bins.split)){
    dfplot              <- as.data.frame(bins.split[[chrom]])
    dfplot$start        <- (-dfplot$start + seqlengths(bins)[chrom])                                                   # Chromosomes start at the top (end is position 1)
    dfplot$counts[which(dfplot$counts >= custom.xlim)] <- custom.xlim                                                  # Change outlier values into quantile 0.999 value
    dfplot.points       <- dfplot[which(dfplot$counts == custom.xlim),]                                                # For plotting outliers (circles)
    ggplt               <- ggplot(dfplot, aes_string(x='start', y='counts'))                                           # ???
    if (!is.null(bins.split[[chrom]]$state)) {                                                                         # If there are states determined
      ggplt             <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts', col='state'), size=0.2)            #   Read count lines
      ggplt             <- ggplt + geom_point(data=dfplot.points,                                                      #   Outliers circles
                             mapping=aes_string(x='start', y='counts', col='state'), size=2, shape=21)
      statelevels       <- unique(c(levels(dfplot$pstate), levels(dfplot$mstate), levels(dfplot$state)))
      ggplt             <- ggplt + scale_color_manual(values=stateColors(statelevels), drop=FALSE)                     #   Colors of the lines and circles corresponding to states. Do not drop levels if not present
    } else {                                                                                                           # If there are no states determined. Is this a possible scenario?
      ggplt             <- ggplt + geom_linerange(aes_string(ymin=0, ymax='counts'), size=0.2, col='gray20')           #   Read count lines
      ggplt             <- ggplt + geom_point(data=dfplot.points,                                                      #   Outliers circles
                            mapping=aes_string(x='start', y='counts'), size=2, shape=21, col='gray20')
    }
    ggplt               <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim,           # Add chromosome backbone as simple rectangle
                             xmin=0, xmax=seqlengths(bins)[chrom], col='white', fill='gray20')
    ggplt               <- ggplt + empty_theme                                                                         # Add theme
    ggplt               <- ggplt + ylab(chrom)                                                                         # Add Chromosome names
    ggplt               <- ggplt + coord_flip(xlim=c(0,max(seqlengths(bins))), ylim=c(-0.6*custom.xlim,custom.xlim))   # Rotate chromosomes 90 degrees. Set x- and y-limits
    ggplts[[chrom]]     <- ggplt
  }
  ncols                 <- ceiling(length(bins.split)/2)                                                               # Combine chromosomes
  plotlist              <- c(rep(list(NULL),ncols), ggplts)
  cowplt                <- plot_grid(plotlist=plotlist, nrow=3, rel_heights=c(2,21,21))                                # Grid with chromosome plots
  cowplt                <- cowplt + cowplot::draw_label(model$ID, x=0.5, y=0.99, vjust=1, hjust=0.5, size=20)          # Add plot title
  return(cowplt)
}

