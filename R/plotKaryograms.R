plotKaryograms <- function(models, cluster=FALSE) {
    if ("aneuHMM" %in% class(models))
        models = list(models)

    segs <- do.call("rbind", lapply(models, function(m) {
        segs <- as.data.frame(m$segments)
        segs$ID <- m$ID
        segs
    }))

    ggplot(segs, aes(color=state), fill="red") +
        geom_segment(aes(x=start, xend=end, y=ID, yend=ID), size=10) +
        facet_grid(. ~ seqnames, scales="free_x") +
        theme(panel.spacing = unit(0, "lines"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
}

