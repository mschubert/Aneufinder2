plotProfile <- function(model) {
    bins <- as.data.frame(model$bins)
    if (!"midpoint" %in% colnames(bins))
        bins$midpoint = (bins$start + bins$end) / 2

    segs <- as.data.frame(model$segments)

    ggplot() +
        geom_point(data=bins, aes(x=midpoint, y=counts)) +
        geom_segment(data=segs, aes(x=start, xend=end, y=mean.counts, yend=mean.counts), size=5) +
        facet_grid(. ~ seqnames, scales="free_x") +
        theme(panel.spacing = unit(0, "lines"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
}
