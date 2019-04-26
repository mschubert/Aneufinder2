stateColors <- function(states=c('zero-inflation', paste0(0:10, '-somy'), 'total')) {
    state.colors        <- c("zero-inflation"="gray90", "0-somy"="gray90","1-somy"="darkorchid3",
                             "2-somy"="springgreen2","3-somy"="red3","4-somy"="gold2",
                             "5-somy"="navy","6-somy"="lemonchiffon","7-somy"="dodgerblue",
                             "8-somy"="chartreuse4","9-somy"="lightcoral","10-somy"="aquamarine2",
                             "total"="black")
    states.with.color   <- intersect(states, names(state.colors))
    cols                <- rep('black', length(states))
    names(cols)         <- states
    cols[states.with.color] <- state.colors[states.with.color]
    return(cols)
}

