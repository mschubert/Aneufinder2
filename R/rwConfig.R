readConfig <- function(configfile=NULL) {
    if (is.null(configfile))
        return(list())

    #TODO: add error handling here
    config <- utils::read.table(configfile, sep="=", fill=TRUE,
                                row.names=NULL, header=FALSE, quote="")

    # Remove spaces
    config[,1] <- gsub(" ","",config[,1])
    config[,2] <- gsub(" ","",config[,2])
    names(config) <- c('argument','value')

    # Remove section headers
    config <- config[config$value != "",]

    configlist <- list() # Turn into config list
    ToParse <- paste0("configlist$", config$argument, " <- ", config$value)

    eval(parse(text=ToParse))
    configlist
}
