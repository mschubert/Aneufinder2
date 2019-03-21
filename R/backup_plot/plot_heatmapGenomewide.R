#' Genome wide heatmap of CNV-state
#'
#' Plot a genome wide heatmap of copy number variation state. This heatmap is
#' best plotted to file, because in most cases it will be too big for cleanly
#' plotting it to screen.
#'
#' @param hmms A list of \code{\link{aneuHMM}} objects or a character vector
#'   with files that contain such objects.
#' @param ylabels A vector with labels for the y-axis. The vector must have the
#'   same length as \code{hmms}. If \code{NULL} the IDs from the
#'   \code{\link{aneuHMM}} objects will be used.
#' @param classes A character vector with the classification of the elements on
#'   the y-axis. The vector must have the same length as \code{hmms}.
#' @param reorder.by.class If \code{TRUE}, the dendrogram will be reordered to
#'   display similar classes next to each other.
#' @param classes.color A (named) vector with colors that are used to
#'   distinguish \code{classes}. Names must correspond to the unique elements in
#'   \code{classes}.
#' @param file A PDF file to which the heatmap will be plotted.
#' @param cluster Either \code{TRUE} or \code{FALSE}, indicating whether the
#'   samples should be clustered by similarity in their CNV-state.
#' @param plot.breakpoints Logical indicating whether breakpoints should be plotted.
#' @param hotspots A \code{\link{GRanges-class}} object with coordinates of
#'   genomic hotspots (see \code{\link{hotspotter}}).
#' @param exclude.regions A \code{\link{GRanges-class}} with regions that will
#'   be excluded from the computation of the clustering. This can be useful to
#'   exclude regions with artifacts.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object or \code{NULL} if a file was specified.
#' @importFrom stats as.dendrogram
#' @importFrom ggdendro dendro_data theme_dendro
#' @importFrom S4Vectors endoapply
#' @export
#' @examples
#' ## Get results from a small-cell-lung-cancer
#' lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#' lung.files <- list.files(lung.folder, full.names=TRUE)
#' ## Get results from the liver metastasis of the same patient
#' liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#' liver.files <- list.files(liver.folder, full.names=TRUE)
#' ## Plot a clustered heatmap
#' classes <- c(rep('lung', length(lung.files)), rep('liver', length(liver.files)))
#' labels <- c(paste('lung',1:length(lung.files)), paste('liver',1:length(liver.files)))
#' heatmapGenomewide(c(lung.files, liver.files), ylabels=labels, classes=classes,
#'                   classes.color=c('blue','red'))
heatmapGenomewide <- function(hmms, ylabels=NULL, classes=NULL,
                              reorder.by.class=TRUE, classes.color=NULL,
                              file=NULL, cluster=TRUE, plot.breakpoints=FALSE,
                              hotspots=NULL, exclude.regions=NULL) {

    ## Check user input
    if (!is.null(ylabels)) {
        if (length(ylabels) != length(hmms)) {
            stop("length(ylabels) must equal length(hmms)")
        }
    }
    if (!is.null(classes)) {
        if (length(classes) != length(hmms)) {
            stop("length(classes) must equal length(hmms)")
        }
    }
    if (length(classes.color)!=length(unique(classes))) {
        stop("'classes.color' must have the same length as unique(classes)")
    }
    if (is.null(names(classes.color))) {
        names(classes.color) <- unique(classes)
    }
    if (!setequal(names(classes.color), unique(classes))) {
        stop("The names of 'classes.color' must be equal to the unique elements in 'classes'")
    }
  if (length(hmms) == 1 & cluster==TRUE) {
    cluster <- FALSE
    warning("Cannot do clustering because only one object was given.")
  }

    ## Load the files
#    hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))

    ## Dataframe with IDs, ylabels and classes
    class.data <- data.frame(ID=sapply(hmms,'[[','ID'))
    class.data$ID <- factor(class.data$ID, levels=class.data$ID)
    if (is.null(ylabels)) {
      class.data$ylabel <- as.character(class.data$ID)
    } else {
      class.data$ylabel <- as.character(ylabels)
    }
    class.data$class <- classes
    
    ## Mapping to match ID with ylabel
    mapping <- class.data$ylabel
    names(mapping) <- class.data$ID

    ## Cluster
    if (reorder.by.class) {
      cl <- clusterHMMs(hmms, cluster=cluster, classes=classes, exclude.regions = exclude.regions)
    } else {
      cl <- clusterHMMs(hmms, cluster=cluster, exclude.regions = exclude.regions)
    }
    hmms <- hmms[cl$IDorder]
    class.data <- class.data[cl$IDorder,]
    class.data$ID <- factor(class.data$ID, levels=class.data$ID)
    ## Extract segements
  segments.list <- GRangesList()
  for (i1 in 1:length(hmms)) {
      hmm <- hmms[[i1]]
      if (is.null(hmm$segments)) {
          segments.list[[hmm$ID]] <- GRanges()
      } else {
            segments.list[[hmm$ID]] <- hmm$segments
      }
  }
  ## Extract breakpoints    
    if (plot.breakpoints) {
      breakpoints <- GRangesList()
      for (i1 in 1:length(hmms)) {
        hmm <- hmms[[i1]]
          if (is.null(hmm$breakpoints)) {
              breakpoints[[hmm$ID]] <- GRanges()
          } else {
              breakpoints[[hmm$ID]] <- hmm$breakpoints
          }
      }
        if (length(breakpoints)==0) {
            plot.breakpoints <- FALSE
        }
    }

    ## Transform coordinates from "chr, start, end" to "genome.start, genome.end"
    ptm <- startTimedMessage("Transforming coordinates ...")
    segments.list <- endoapply(segments.list, transCoord)
    if (plot.breakpoints) {
        breakpoints <- endoapply(breakpoints, transCoord)
    }
    stopTimedMessage(ptm)

    ## Data.frame for plotting
    ptm <- startTimedMessage("Making the plot ...")
    df <- list()
    for (i1 in 1:length(segments.list)) {
        df[[length(df)+1]] <- data.frame(start=segments.list[[i1]]$start.genome, end=segments.list[[i1]]$end.genome, seqnames=seqnames(segments.list[[i1]]), ID=names(segments.list)[i1], state=segments.list[[i1]]$state)
    }
    df <- do.call(rbind, df)
    df$ID <- factor(df$ID, levels=levels(class.data$ID))
    df$ylabel <- mapping[as.character(df$ID)]
    
    if (plot.breakpoints) {
        df.breakpoints <- list()
        for (i1 in 1:length(breakpoints)) {
          if (length(breakpoints[[i1]]) > 0) {
                df.breakpoints[[length(df.breakpoints)+1]] <- data.frame(start=breakpoints[[i1]]$start.genome, end=breakpoints[[i1]]$end.genome, seqnames=seqnames(breakpoints[[i1]]), ID=names(segments.list)[i1], mid=(breakpoints[[i1]]$start.genome + breakpoints[[i1]]$end.genome)/2)
          } else {
                df.breakpoints[[length(df.breakpoints)+1]] <- data.frame(start=numeric(), end=numeric(), seqnames=character(), ID=character(), mid=numeric())
          }
        }
        df.breakpoints <- do.call(rbind, df.breakpoints)
      df.breakpoints$ID <- factor(df.breakpoints$ID, levels=levels(class.data$ID))
      df.breakpoints$ylabel <- mapping[as.character(df.breakpoints$ID)]
    }
    
    # Chromosome lines
    cum.seqlengths <- cumsum(as.numeric(seqlengths(segments.list[[1]])))
    names(cum.seqlengths) <- seqlevels(segments.list[[1]])
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(segments.list[[1]])
    label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(segments.list[[1]]) )
    df.chroms <- data.frame(y=c(0,cum.seqlengths), x=1, xend=length(segments.list))

    ### Plot ###
    pltlist <- list()
    widths <- vector()

    ## Prepare the plot
    df$state <- factor(df$state, levels=names(sort(initializeStates(levels(df$state))$multiplicity)))
    df$x <- as.numeric(df$ID) # transform all x-coordiantes to numeric because factors and numerics get selected different margins
    ggplt <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='x', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + scale_x_continuous(name="sample", breaks=1:length(unique(df$ylabel)), labels=unique(df$ylabel))
    ggplt <- ggplt + scale_color_manual(values=stateColors(levels(df$state)))
    ggplt <- ggplt + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20), axis.line=element_blank(), axis.title.x=element_blank())
    # ggplt <- ggplt + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
    ggplt <- ggplt + geom_segment(aes_string(x='x', xend='xend', y='y', yend='y'), data=df.chroms, col='black')
    ggplt <- ggplt + coord_flip()
    if (plot.breakpoints) {
        df.breakpoints$x <- as.numeric(df.breakpoints$ID)
        ggplt <- ggplt + geom_linerange(data=df.breakpoints, mapping=aes_string(x='x', ymin='start', ymax='end'), size=2) + ylab('') + geom_point(data=df.breakpoints, mapping=aes_string(x='x', y='mid'))
    }
    if (!is.null(hotspots)) {
      if (length(hotspots) > 0) {
          df.hot <- as.data.frame(transCoord(hotspots))
          df.hot$xmin <- 0
          df.hot$xmax <- length(class.data$ID)+1
          ggplt <- ggplt + geom_rect(data=df.hot, mapping=aes_string(xmin='xmin', xmax='xmax', ymin='start.genome', ymax='end.genome', alpha='num.events'), fill='hotpink4') + scale_alpha_continuous(name='breakpoints', range=c(0.4,0.8))
      }
    }
    width.heatmap <- sum(as.numeric(seqlengths(hmms[[1]]$bins))) / 3e9 * 150 # human genome (3e9) roughly corresponds to 150cm
    height <- max(length(hmms) * 0.5, 2)
    pltlist[['heatmap']] <- ggplt
    widths['heatmap'] <- width.heatmap
    ## Make the classification bar
    if (!is.null(classes)) {
        width.classes <- 5
      class.data$x <- as.numeric(class.data$ID)  # transform all x-coordiantes to numeric because factors and numerics get selected different margins
        ggclass <- ggplot(class.data) + geom_linerange(aes_string(ymin=0, ymax=1, x='x', col='class'), size=5) + guides(col=FALSE) + xlab("class")
      ggclass <- ggclass + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title.x=element_blank())
        ggclass <- ggclass + coord_flip()
        if (!is.null(classes.color)) {
            ggclass <- ggclass + scale_color_manual(breaks=names(classes.color), values=classes.color)
        }
        pltlist[['classbar']] <- ggclass
        widths['classbar'] <- width.classes
    }
    ## Prepare the dendrogram
    if (!is.null(cl$hclust)) {
        dhc <- stats::as.dendrogram(cl$hclust)
        ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
        ggdndr <- ggplot(ddata$segments) + geom_segment(aes_string(x='x', xend='xend', y='y', yend='yend')) + scale_y_reverse()
        ggdndr <- ggdndr + coord_flip()
        ggdndr <- ggdndr + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title=element_blank())
        width.dendro <- 20
        pltlist[['dendro']] <- ggdndr
        widths['dendro'] <- width.dendro
    }
    cowplt <- cowplot::plot_grid(plotlist=rev(pltlist), align='h', ncol=length(pltlist), rel_widths=rev(widths))
    stopTimedMessage(ptm)

    ## Plot to file
    if (!is.null(file)) {
        ptm <- startTimedMessage("Plotting to file ",file," ...")
        ggsave(file, cowplt, width=sum(widths), height=height, units='cm', limitsize=FALSE)
        stopTimedMessage(ptm)
    } else {
        return(cowplt)
    }

}
