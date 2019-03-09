checkClass <- function(conf) {
    # [FOLDERS]
    if (!is.character(conf[['inputfolder']]))
        stop("Argument isn't a character: 'inputfolder'")
    if (!is.character(conf[['outputfolder']]))
        stop("Argument isn't a character: 'outputfolder'")

    # [GENERAL]
    if (!is.numeric(conf[['numCPU']]))
        stop("Argument isn't numeric: 'numCPU'")
    if (!is.logical(conf[['reuse.existing.files']]))
        stop("Argument isn't logical (TRUE or FALSE): 'reuse.existing.files'")

    # [BINNING]
    if(!is.numeric(conf[['binsizes']]))
        stop("Argument isn't numeric: 'binsizes'")
    if (!is.numeric(conf[['stepsizes']]) && !is.null(conf[['stepsizes']]))
        stop("Argument isn't numeric or null: 'stepsizes'")
    if (!is.character(conf[['variable.width.reference']]) &&
            !is.null(conf[['variable.width.reference']]))
        stop("Argument isn't a character or null: 'variable.width.reference'")
    if (!is.numeric(conf[['reads.per.bin']]) && !is.null(conf[['reads.per.bin']]))
        stop("Argument isn't numeric or null: 'reads.per.bin'")
    if (!is.logical(conf[['pairedEndReads']]))
        stop("Argument isn't logical (TRUE or FALSE): 'pairedEndReads'")
    if (!is.character(conf[['assembly']]) && is.data.frame(conf[['assembly']]) &&
            !is.null(conf[['assembly']]))
        stop("Argument isn't a character, a data.frame or null: 'assembly'")
    if (!is.numeric(conf[['chromosomes']]) && !is.character(conf[['chromosomes']]))
        stop("Argument isn't numeric or a character: 'chromosomes'")
    if (!is.logical(conf[['remove.duplicate.reads']]))
        stop("Argument isn't logical (TRUE or FALSE): 'remove.duplicate.reads'")
    if (!is.numeric(conf[['min.mapq']]))
        stop("Argument isn't numeric: 'min.mapq'")
    if (!is.character(conf[['blacklist']]) && (class(blacklist) != 'GRanges') &&
            !is.null(conf[['blacklist']]))
        stop("Argument isn't a character, 'GRanges' or null: 'blacklist'")
    if (!is.logical(conf[['use.bamsignals']]))
        stop("Argument isn't logical (TRUE or FALSE): 'use.bamsignals'")
    if (!is.logical(conf[['reads.store']]))
        stop("Argument isn't logical (TRUE or FALSE): 'reads.store'")
  
    # [CORRECTION]
    if (!is.character(conf[['correction.method']]) &&
            !is.null(conf[['correction.method']]))
        stop("Argument isn't a character or null: 'correction.method'")
    if (!is.character(conf[['GC.BSgenome']]) &&
            (class(conf[['GC.BSgenome']]) != 'BSgenome') &&
            !is.null(conf[['GC.BSgenome']]))
        stop("Argument isn't a character, 'BSgenome' or null: 'GC.BSgenome'")

    # [COPYNUMBERCALLING]
    if (!is.character(conf[['method']]))
        stop("Argument isn't a character: 'method'")
    if (!is.logical(conf[['strandseq']]))
        stop("Argument isn't logical (TRUE or FALSE): 'strandseq'")
    if (!is.numeric(conf[['R']]))
        stop("Argument isn't numeric: 'R'")
    if (!is.numeric(conf[['sig.lvl']]))
        stop("Argument isn't numeric: 'sig.lvl'")
    if (!is.numeric(conf[['eps']]))
        stop("Argument isn't numeric: 'eps'")
    if (!is.numeric(conf[['max.time']]))
        stop("Argument isn't numeric: 'max.time'")
    if (!is.numeric(conf[['max.iter']]))
        stop("Argument isn't numeric: 'max.iter'")
    if (!is.numeric(conf[['num.trials']]))
        stop("Argument isn't numeric: 'num.trials'")
    if (!is.character(conf[['states']]))
        stop("Argument isn't a character: 'states'")
    if (!is.character(conf[['most.frequent.state']]))
        stop("Argument isn't a character: 'most.frequent.state'")
    if (!is.character(conf[['most.frequent.state.strandseq']]))
        stop("Argument isn't a character: 'most.frequent.state.strandseq'")
}
