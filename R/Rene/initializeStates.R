
# ---------------------------------------------------------------------------------------------------------------------
#                                                  initializeStates
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 09-11-18
# Last modified: 09-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

initializeStates <- function(states){
  somy.states <- grep('somy', states, value=TRUE)
  somy.numbers <- as.integer(sapply(strsplit(somy.states, '-somy'), '[[', 1))
  names(somy.numbers) <- somy.states
  
  if("zero-inflation" %in% states){
    multiplicity <- c("zero-inflation"=0, somy.numbers)
  }else{
    multiplicity <- somy.numbers
  }
  
  levels.distributions <- c('delta','dgeom','dnbinom','dbinom')
  distributions <- rep(NA, length(states))
  names(distributions) <- states
  distributions[states=='zero-inflation'] <- 'delta'
  distributions[states=='0-somy'] <- 'dgeom'
  distributions[(states != 'zero-inflation') & (states != '0-somy')] <- 'dnbinom'
  states <- factor(states, levels=states)
  distributions <- factor(distributions, levels=levels.distributions)
  l <- list(states=states, distributions=distributions, multiplicity=multiplicity)
  return(l)
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------