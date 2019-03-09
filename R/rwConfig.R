
# ---------------------------------------------------------------------------------------------------------------------
#                                                      rwConfig
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 01-11-18
# Last modified: 01-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_readConfig <- function(configfile){
  config        <- utils::read.table(configfile, sep="=", fill=TRUE, row.names=NULL, header=FALSE, quote="")
  config[,1]    <- gsub(" ","",config[,1])                                                                             # Remove spaces
  config[,2]    <- gsub(" ","",config[,2])
  names(config) <- c('argument','value')
  config        <- config[config$value != "",]                                                                         # Remove section headers
  configlist    <- list()                                                                                              # Turn into config list
  ToParse       <- paste0("configlist$", config$argument, " <- ", config$value)
  eval(parse(text=ToParse)) 
  return(configlist) 
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------