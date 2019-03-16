
# =====================================================================================================================
#                                                    Aneufinder.R
# =====================================================================================================================
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 25-09-18
# Last modified: 07-03-19
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------


# =====================================================================================================================
# PART 0 of 4 | ARGUMENTS AND DEFAULT SETTINGS
# =====================================================================================================================

Aneufinder <- function(inputfolder,
                       outputfolder,
                       configfile=NULL,
                       numCPU=1,                                                                                       # [GENERAL]
                       reuse.existing.files=TRUE,
                       binsizes=1e6,                                                                                   # [BINNING]
                       stepsizes=binsizes,
                       variable.width.reference=NULL,
                       reads.per.bin=NULL,
                       pairedEndReads=FALSE,
                       assembly=NULL,
                       chromosomes=NULL,
                       remove.duplicate.reads=TRUE,
                       min.mapq=10,
                       blacklist=NULL,
                       use.bamsignals=FALSE,
                       reads.store=FALSE,
                       correction.method=NULL,                                                                         # [CORRECTION]
                       GC.BSgenome=NULL,
                       method='edivisive',                                                                             # [COPYNUMBERCALLING]
                       strandseq=FALSE,
                       R=10,
                       sig.lvl=0.1,
                       eps=0.01,
                       max.time=60,
                       max.iter=5000,
                       num.trials=15,
                       states=c('zero-inflation',paste0(0:10,'-somy')),
                       most.frequent.state='2-somy',                                                                   # N: New!
                       most.frequent.state.strandseq='1-somy',                                                         # N: New!
                       confint=NULL,
                       refine.breakpoints=FALSE,
                       hotspot.bandwidth=NULL,
                       hotspot.pval=5e-2,
                       cluster.plots=TRUE){


# =====================================================================================================================
# PART 1 OF 4 | PREPARATION
# =====================================================================================================================

# ---------------------------------------------------------------------------------------------------------------------
# PART 1.01 | READ CONFIG FILE
# ---------------------------------------------------------------------------------------------------------------------

conf                    <- NULL

if(!is.null(configfile)){
  errstring <- tryCatch({
    conf                <- RW_readConfig(configfile)
    errstring           <- ''
  }, error = function(err) {
    errstring           <- paste0("Could not read configuration file ",configfile)
  })
  if(errstring!='') {
    stop(errstring)
  }
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.02 | COMBINE ARGUMENTS CONFIG FILE WITH SUPPLIED ARGUMENTS
# ---------------------------------------------------------------------------------------------------------------------

params                   <- list(inputfolder=inputfolder, outputfolder=outputfolder, numCPU=numCPU,
                                 reuse.existing.files=reuse.existing.files, binsizes=binsizes, stepsizes=stepsizes,
                                 variable.width.reference=variable.width.reference, reads.per.bin=reads.per.bin,
                                 pairedEndReads=pairedEndReads, assembly=assembly, chromosomes=chromosomes,
                                 remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, blacklist=blacklist,
                                 reads.store=reads.store, use.bamsignals=use.bamsignals,
                                 correction.method=correction.method, GC.BSgenome=GC.BSgenome, method=method,
                                 strandseq=strandseq, eps=eps, max.time=max.time, max.iter=max.iter,
                                 num.trials=num.trials, states=states, most.frequent.state=most.frequent.state,
                                 most.frequent.state.strandseq=most.frequent.state.strandseq, R=R, sig.lvl=sig.lvl,
                                 confint=confint, refine.breakpoints=refine.breakpoints,
                                 hotspot.bandwidth=hotspot.bandwidth, hotspot.pval=hotspot.pval,
                                 cluster.plots=cluster.plots)
conf                    <- c(conf, params[setdiff(names(params),names(conf))])


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.03 | CHECK CLASS
# ---------------------------------------------------------------------------------------------------------------------

checkClass(conf=conf)


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.04 | GET INPUT FILES AND CHECK FORMAT
# ---------------------------------------------------------------------------------------------------------------------

input.files             <- list.files(inputfolder, full.names=TRUE, pattern='\\.bam$|\\.bed$|\\.bed\\.gz$')

if(length(input.files) == 0){
  stop("None of the input files have the correct format. Expected formats are '.bam', '.bed' and '.bed.gz'")
}

files.clean             <- sub('\\.gz$','',input.files)
input.type              <- unique(sapply(strsplit(files.clean,'\\.'),function(x){rev(x)[1]}))

if(length(input.type) == 2){
  stop("Both bam and bed files in input directory. Only one type allowed.")                                            # N: Check if input files is a mix of different types (bam, bed, bed.gz). Allow only one.
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.05 | [GENERAL]
# ---------------------------------------------------------------------------------------------------------------------

if(conf[['reuse.existing.files']] == FALSE){                                                                           # N: Delete old directory if desired
  if(file.exists(outputfolder)){
    message("Deleting old directory: ",outputfolder)
    unlink(outputfolder, recursive=TRUE)
  }
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.06 | [BINNING]
# ---------------------------------------------------------------------------------------------------------------------

if(is.null(conf[['stepsizes']])){
  conf[['stepsizes']]   <- conf[['binsizes']]
}

if(length(conf[['binsizes']]) != length(conf[['stepsizes']])){
  stop("Need one element in 'stepsizes' for each element in 'binsizes'.")
}

if(any(conf[['binsizes']] < conf[['stepsizes']])){
  stop("'stepsizes' must be smaller/equal than 'binsizes'")
}

if(!is.null(conf[['variable.width.reference']])){
  file.clean            <- sub('\\.gz$','',conf[['variable.width.reference']])
  ref.type              <- rev(strsplit(file.clean,'\\.')[[1]])[1]
  if((ref.type != 'bam') & (ref.type != 'bed')){
    stop("The variable width reference file does not have the correct format.
    The expected formats are '.bam', '.bed' and '.bed.gz'")
  }
  if(!file.exists(conf[['variable.width.reference']])){
    stop("variable.width.reference file '",conf[['variable.width.reference']],"' does not exist.")
  }
}

if(!is.null(conf[['reads.per.bin']])){
  if(conf[['reads.per.bin']] < 1){
    stop("The number of reads per bin is smaller than the minimum allowed (<1).")
  }
}

if(conf[['min.mapq']] < 0){
  stop("Unusual low 'min.mapq': ",conf[['min.mapq']]) 
}

if(is.character(conf[['blacklist']])){
  file.clean            <- sub('\\.gz$','',conf[['blacklist']])
  black.type            <- rev(strsplit(file.clean,'\\.')[[1]])[1]
  if(black.type != 'bed'){
    stop("The blacklist has the wrong file format: ",file.format,". Allowed formats are: 'bed' (tab delimited)")
  }
  if(!file.exists(conf[['blacklist']])){
    stop("Blacklist file '",conf[['blacklist']],"' does not exist.")
  }
  conf[['blacklist']]   <- suppressMessages(readBed(conf[['blacklist']], track.line="auto", remove.unusual=FALSE,
                             zero.based=TRUE))
  if(class(conf[['blacklist']])[1] != "GRanges"){
    stop("Something went wrong with reading the bed file. Please check whether the bed file is tab delimited.")
  }
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.07 [CORRECTION]
# ---------------------------------------------------------------------------------------------------------------------

if(!all(conf[['correction.method']] %in% c('GC'))){                                                                    # Check if the correct correction method has been given.
  stop("Unknown correction method: ",paste(setdiff(conf[['correction.method']],c("GC")),collapse=', '),".
    Allowed methods are: 'GC'.")
}

if('GC' %in% conf[['correction.method']] & is.null(conf[['GC.BSgenome']])){                                            # Check whether method 'GC' is given in combination with 'GC.genome'
  stop("Option 'GC.bsgenome' has to be given if correction.method='GC'.")
}

if(!is.null(conf[['GC.BSgenome']])){
  if(is.character(conf[['GC.BSgenome']])){
    eval(parse(text=paste0("conf[['GC.BSgenome']] <- ",conf[['GC.BSgenome']])))
  }else if(class(conf[['GC.BSgenome']]) != 'BSgenome'){
    stop("Unknown class for 'GC.BSgenome': ",class(GC.BSgenome),".
      'GC.BSgenome' should either be of class 'character' or 'BSgenome'")
  }
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.08 [COPYNUMBERCALLING]
# ---------------------------------------------------------------------------------------------------------------------

if(!all(conf[['method']] %in% c('HMM','dnacopy','edivisive'))){ 
  stop("Unknown copynumber calling method ('method'): ",paste(setdiff(conf[['method']],c('HMM','dnacopy','edivisive')),
    collapse=', '),". Allowed methods are: 'HMM', 'dnacopy' and 'edivisive'.")
}

if(conf[['R']] < 1){
  stop("The maximum number of random permutations to use in each iteration of the permutation test should be 1 or
    higher (method: edivisive). Current value ('R'): ",conf[['R']])
}

if(conf[['sig.lvl']] <= 0 | conf[['sig.lvl']] > 1){
  stop("The statistical significance level for a proposed change point should be between 0 and 1. 
   Current value ('sig.lvl'): ",conf[['sig.lvl']])
}

if(conf[['eps']] <= 0){
  stop("The Convergence threshold for the Baum-Welch algorithm should be higher than 0.
    Current value ('eps'): ",conf[['eps']])
}

if(conf[['eps']] > 0.1){
  warning("Unusual high number for the Convergence threshold for the Baum-Welch algorithm (>0.1).
    Current value ('eps'): ",conf[['eps']])
}

conf[['max.time']]      <- round(conf[['max.time']])                                                                   # Ensures that we get an integer.

if(conf[['max.time']] < -1 | conf[['max.time']] == 0){
  stop("The maximum running time in seconds for the Baum-Welch algorithm can have a value higher than zero or -1
    (no limit). Current value ('max.time'): ",conf[['max.time']])
}

conf[['max.iter']]      <- round(conf[['max.iter']])                                                                   # Ensures that we get an integer.

if(conf[['max.iter']] < -1 | conf[['max.iter']] == 0){
  stop("The maximum number of iterations for the Baum-Welch algorithm can have a value higher than zero or -1
    (no limit). Current value ('max.iter'): ",conf[['max.iter']])
}

conf[['num.trials']]    <- round(conf[['num.trials']])                                                                 # Ensures that we get an integer.

if(conf[['num.trials']] <= 0){
  stop("The number of trials to find a fit where state 'most.frequent.state' is most frequent should be higher than 0.
    Current value ('num.trials'): ",conf[['num.trials']])
}

if(any(!grepl('zero-inflation|^[0-9]+-somy|+[0-9]+-somy',conf[['states']]))){
  stop("")
}
if(any(table(conf[['states']]) != 1)){                                                                                 # Check non-unique states.
  stop("States are not unique.")
}
if(any(grepl('zero-inflation',conf[['states']]))){
  if(grep('zero-inflation',conf[['states']]) != 1){
    stop("The zero-inflation state should be the first of all states.")
  }
}
state.somy              <- grep('-somy',conf[['states']],value=TRUE)
state.num               <- substr(state.somy,1,nchar(state.somy)-5)
state.plus              <- grep('^\\+',state.num)
if(length(state.plus) > 0){
  if(length(state.plus) > 1){
    stop("There is more than one +[number]-somy state.")
  }
  if(state.plus == length(state.num)){
    state.num           <- state.num[-state.plus]
  }
}
if(any(state.num != sort(as.numeric(state.num)))){
  stop("States are not ordered.")
}

if(!conf[['most.frequent.state']] %in% conf[['states']]){
  stop("argument 'most.frequent.state' must be one of c(",paste(states, collapse=","),")")
}

if(!conf[['most.frequent.state.strandseq']] %in% conf[['states']]){
  stop("argument 'most.frequent.state.strandseq' must be one of c(",paste(states, collapse=","),")")
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.09 | CHECK CHROMOSOME FORMATS AND OVERLAP WITH SPECIFIED CHROMOSOMES
# ---------------------------------------------------------------------------------------------------------------------

if(input.type == 'bam'){                                                                                               # N: We take the chromosome format of the input files as the format to work with.
  seq.names             <- GenomeInfoDb::seqlevels(Rsamtools::BamFile(input.files[1]))
}else if(input.type == 'bed' | input.type == 'bed.gz'){
  seq.names             <- GenomeInfoDb::seqlevels(readBed(input.files[1]))
}

if(all(grepl('^chr',seq.names))){
  chrom.type            <- 'chr'
}else if(all(!grepl('^chr',seq.names))){
  chrom.type            <- 'num'
}else{                                                                                                                 # Q: Check if there is a mix of chr and not chr?
  stop("Inconsistency in chromosome names input files. Some start with 'chr' while others do not.")
}

chrom.missing           <- setdiff(conf[['chromosomes']],seq.names)

if(length(chrom.missing) == length(conf[['chromosomes']])){
  chr.string            <- paste0(chrom.missing, collapse=', ')
  stop("The specified chromosomes ",chr.string, " are not found in the data (sam/bed files).
    Pay attention to the naming convention in your data, e.g. 'chr1' or '1'.")
}else if(length(chrom.missing) > 0){
  chr.string            <- paste0(chrom.missing, collapse=', ')
  warning(paste0('Not using chromosomes ',chr.string,' because they are not found in the data (sam/bed files).'))
}

if(!is.null(conf[['variable.width.reference']])){
  if(ref.type == 'bam'){                                                                                               # N: Check 'variable.width.reference'
    seq.names           <- GenomeInfoDb::seqlevels(Rsamtools::BamFile(conf[['variable.width.reference']]))
  }else if(ref.type == 'bed'){
    seq.names           <- GenomeInfoDb::seqlevels(readBed(conf[['variable.width.reference']]))
  }
  if(all(grepl('^chr',seq.names))){
    if(chrom.type == 'num'){
      stop("The specified chromosomes do not exist in the data.
        Pay attention to the naming convention in your data, e.g. 'chr1' or '1'.")
    }
  }else if(all(!grepl('^chr',seq.names))){
    if(chrom.type == 'chr'){
      stop("The specified chromosomes do not exist in the data.
        Pay attention to the naming convention in your data, e.g. 'chr1' or '1'.")
    }
  }else{                                                                                                               # Q: Check if there is a mix of chr and not chr?
    stop("Inconsistency in chromosome names variable width reference. Some start with 'chr' while others do not.")
  }
  chrom.missing         <- setdiff(conf[['chromosomes']],seq.names)
  if(length(chrom.missing) == length(conf[['chromosomes']])){
    chr.string          <- paste0(chrom.missing, collapse=', ')
    stop("The specified chromosomes ",chr.string, " are not found in the data (variable.width.reference).
      Pay attention to the naming convention in your data, e.g. 'chr1' or '1'.")
  }else if(length(chrom.missing) > 0){
    chr.string          <- paste0(chrom.missing, collapse=', ')
    warning(paste0('Not using chromosomes ',chr.string,' because they are not found in the data
     (variable.width.reference).'))
  }
}

if(!is.null(conf[['blacklist']])){
  seq.names             <- GenomeInfoDb::seqlevels(conf[['blacklist']])
  if(all(grepl('^chr',seq.names))){
    if(chrom.type == 'num'){
      seqlevels(conf[['blacklist']]) <- sub('chr','',seqlevels(conf[['blacklist']]))
    }
  }else if(all(!grepl('^chr',seq.names))){
    if(chrom.type == 'chr'){
      seqlevels(conf[['blacklist']]) <- paste0('chr',seqlevels(conf[['blacklist']]))
    }
  }else{                                                                                                               # Q: Check if there is a mix of chr and not chr?
    stop("Inconsistency in chromosome names blacklist. Some start with 'chr' while others do not.")
  }
  chrom.missing         <- setdiff(conf[['chromosomes']],seq.names)
  if(length(chrom.missing) == length(conf[['chromosomes']])){
    chr.string          <- paste0(chrom.missing, collapse=', ')
    stop("The specified chromosomes ",chr.string," are not found in the data (blacklist).
      Pay attention to the naming convention in your data, e.g. 'chr1' or '1'.")
  }else if(length(chrom.missing) > 0){
    chr.string          <- paste0(chrom.missing, collapse=', ')
    warning(paste0("Not using chromosomes ",chr.string," because they are not found in the data (blacklist)."))
  }
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.10 | GET CHROMOSOME LENGTHS
# ---------------------------------------------------------------------------------------------------------------------

chrom.lengths           <- RW_chromosomeLengths(assembly=conf[['assembly']],chrom.type=chrom.type,                     # N: We can indicate what kind of source is used.
                             input.files=input.files)
chrom.lengths           <- chrom.lengths[which(names(chrom.lengths) %in% conf[['chromosomes']])]
chrom.missing           <- setdiff(conf[['chromosomes']],names(chrom.lengths))

if(length(chrom.missing) == length(conf[['chromosomes']])){
  chr.string            <- paste0(chrom.missing, collapse=', ')
  stop("The specified chromosomes ",chr.string," are not found within the object or file that is specifying the
    chromosome lengths. Pay attention to the naming convention in your data, e.g. 'chr1' or '1'.")
}else if(length(chrom.missing) > 0){
  chr.string            <- paste0(chrom.missing, collapse=', ')
  warning(paste0('Not using chromosomes ',chr.string,' because the object or file that is specifying the chromosome
    lengths does not contain them.'))
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 1.11 | CREATE OUTPUT DIRECTORY
# ---------------------------------------------------------------------------------------------------------------------

if(!file.exists(conf[['outputfolder']])){
  dir.create(conf[['outputfolder']])
}


# =====================================================================================================================
# PART 2 OF 4 | FILTERING, BINNING AND CORRECTING THE DATA
# =====================================================================================================================

# ---------------------------------------------------------------------------------------------------------------------
# PART 2.01 | MAKE BIN LIST
# ---------------------------------------------------------------------------------------------------------------------

bins.list               <- RW_fixedWidthBins(chrom.lengths=chrom.lengths, binsizes=conf[['binsizes']],
                             stepsizes=conf[['stepsizes']])

if(!is.null(conf[['variable.width.reference']])){                                                                      # Q: This part takes a long time. Can we skip it if the binned files are already present? Check if all expected files are there based on names input files etc.
  if(ref.type == 'bam'){
    bam.index.ref       <- paste0(conf[['variable.width.reference']],".bai")
    if(!file.exists(bam.index.ref)){
      bam.index.ref     <- Rsamtools::indexBam(conf[['variable.width.reference']])
      warning("Couldn't find BAM index-file. Creating our own file ", bam.index.ref," instead.")
    }
    reads.ref           <- RW_bam2GRanges(bamfile=conf[['variable.width.reference']], bamindex=bam.index.ref,
                             chrom.lengths=chrom.lengths, pairedEndReads=conf[['pairedEndReads']],
                             remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']],
                             blacklist=blacklist)
  }else if(ref.type == 'bed'){
    reads.ref           <- RW_bed2GRanges(bedfile=conf[['variable.width.reference']], chrom.lengths=chrom.lengths,
                             remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']],
                             blacklist=blacklist)
  }
  binned.ref            <- list()
  for(cbss in names(bins.list)){                                                                                       # N: cbss: combination bin step and size
    binned.ref[[cbss]]  <- RW_binReads(reads=reads.ref, bins=bins.list[[cbss]])
  }
  #binned.ref            <- RW_binReads(reads=reads.ref, bins.list=bins.list)
  
  bins.list             <- NULL
  bins.list             <- RW_variableWidthBins(reads=reads.ref, binned.list=binned.ref)
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 2.02 | GET READ DATA AND FILTER
# ---------------------------------------------------------------------------------------------------------------------

path.filtered.reads     <- file.path(conf[['outputfolder']],'filtered')
if(!file.exists(path.filtered.reads)){dir.create(path.filtered.reads)}

files.to.do             <- setdiff(basename(input.files),list.files(path.filtered.reads,full.names=FALSE))
files.to.do             <- file.path(conf[['inputfolder']],files.to.do)

for(file.cur in files.to.do){
  if(input.type == "bam"){
    bam.index           <- paste0(file.cur,".bai")
    if(!file.exists(bam.index)){
      bam.index         <- Rsamtools::indexBam(file.cur)
      warning("Couldn't find BAM index-file. Creating our own file ",bam.index," instead.")
    }
    reads               <- RW_bam2GRanges(bamfile=file.cur, bamindex=bam.index, chrom.lengths=chrom.lengths,
                             pairedEndReads=conf[['pairedEndReads']],
                             remove.duplicate.reads=conf[['remove.duplicate.reads']],
                             min.mapq=conf[['min.mapq']], blacklist=blacklist)
  }else if(input.type == "bed"){
    reads               <- RW_bed2GRanges(bedfile=file.cur, chrom.lengths=chrom.lengths,
                             remove.duplicate.reads=conf[['remove.duplicate.reads']],
                             min.mapq=conf[['min.mapq']], blacklist=blacklist)
  }
  save(reads,file=file.path(path.filtered.reads,paste0(basename(file.cur),'.Rdata')))
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 2.03 | BIN THE READS
# ---------------------------------------------------------------------------------------------------------------------

path.uncorrected.bins   <- file.path(conf[['outputfolder']],'binned')
if(!file.exists(path.uncorrected.bins)){dir.create(path.uncorrected.bins)}

files.to.do             <- list.files(path.filtered.reads,full.names=TRUE)

for(file.cur in files.to.do){
  reads                 <- get(load(file.cur))
  for(ibss in 1:length(conf[['binsizes']])){                                                                           # ibss: index bin step size combinations
    binsize             <- conf[['binsizes']][ibss]
    stepsize            <- conf[['stepsizes']][ibss]
    combi               <- paste0("binsize_",format(binsize,scientific=TRUE,trim=TRUE),"_stepsize_",
                             format(stepsize,scientific=TRUE,trim=TRUE))
    inp_file            <- basename(file.cur)
    inp_file            <- substr(inp_file,1,(nchar(inp_file)-6))
    file.save           <- file.path(path.uncorrected.bins,paste0(inp_file,"_",combi,".RData"))
    if(!file.exists(file.save)){
      binned            <- RW_binReads(reads=reads,bins=bins.list[[combi]])
      save(binned,file=file.save)
    }
  }
}


# ---------------------------------------------------------------------------------------------------------------------
# PART 2.04 | CORRECT READ COUNT
# ---------------------------------------------------------------------------------------------------------------------

if(!is.null(conf[['correction.method']])){                                                                             # Q: Note: after correction the total bin count differs from the sum of positive and negative strand due to rounding. Should we somehow 'fix' this?
  path.corrected.bins   <- paste0(path.uncorrected.bins,'-',conf[['correction.method']])
  if(!file.exists(path.corrected.bins)){dir.create(path.corrected.bins)}
  if(conf[['correction.method']] == 'GC'){
    bins.list.GC        <- RW_getGCContentBins(bins.list=bins.list,GC.BSgenome=conf[['GC.BSgenome']])
  }
  files.to.do           <- setdiff(list.files(path.uncorrected.bins,full.names=FALSE),
                             list.files(path.corrected.bins,full.names=FALSE))
  files.to.do           <- file.path(path.uncorrected.bins,files.to.do)
  for(file.cur in files.to.do){
    if(conf[['correction.method']] == 'GC'){
      binned            <- get(load(file.cur))
      split_res         <- strsplit(basename(file.cur),"_binsize_")                                                    # N: Split on "_binsize_"
      combi             <- paste0("binsize_",substr(split_res[[1]][2],1,(nchar(split_res[[1]][2])-6)))                 # N: Remove ".RData" and add again "binsize_"
      binned.GC         <- merge(binned,bins.list.GC[[combi]])
      binned.GC.cor     <- RW_correctGC(binned.gc=binned.GC,method='loess')
      save(binned.GC.cor,file=file.path(path.corrected.bins,basename(file.cur)))
    }
  }
}else{
  path.corrected.bins   <- path.uncorrected.bins  
}


# =====================================================================================================================
# PART 3 OF 4 | RUN MODELS
# =====================================================================================================================

# ---------------------------------------------------------------------------------------------------------------------
# PART 3.01 | RUN MODELS
# ---------------------------------------------------------------------------------------------------------------------

path.model              <- file.path(conf[['outputfolder']],'MODELS')
if(!file.exists(path.model)){dir.create(path.model)}

####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE
####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE
####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE

readspath               <- file.path(outputfolder,'data')
modelpath               <- file.path(outputfolder, 'MODELS')
refinedmodelpath        <- paste0(modelpath, '_refined')
plotpath                <- file.path(outputfolder, 'PLOTS')
browserpath             <- file.path(outputfolder, 'BROWSERFILES')


numcpu                  <- conf[['numCPU']]


binsizes                <- conf[['binsizes']]
stepsizes               <- conf[['stepsizes']]
reads.per.bins          <- conf[['reads.per.bin']]

patterns                <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=TRUE, trim=TRUE),'_stepsize_',format(stepsizes, scientific=TRUE, trim=TRUE),'_'))
patterns                <- setdiff(patterns, c('reads.per.bin__','binsize__'))

patterns                <- substr(patterns,1,(nchar(patterns)-1))


cl <- parallel::makeCluster(numcpu)
doParallel::registerDoParallel(cl)
on.exit(
  if (conf[['numCPU']] > 1) {
    parallel::stopCluster(cl)
  }
)



#############################################################################################################################


for(method in conf[['method']]){
  path.method           <- file.path(path.model,paste0('method-',method))
  if(!file.exists(path.method)){dir.create(path.method)}
  files.to.do           <- setdiff(list.files(path.corrected.bins,full.names=FALSE),
                             list.files(path.method,full.names=FALSE))
  files.to.do           <- file.path(path.corrected.bins,files.to.do)
  for(file.cur in files.to.do){
    binned              <- get(load(file.cur))
    model               <- RW_findCNVs(strandseq=conf[['strandseq']], binned=binned, ID=attr(binned,'ID'), method=method,
                             R=conf[['R']], sig.lvl=conf[['sig.lvl']], eps=conf[['eps']], max.time=conf[['max.time']],
                             max.iter=conf[['max.iter']], num.trials=conf[['num.trials']], states=conf[['states']],
                             most.frequent.state=conf[['most.frequent.state']],
                             most.frequent.state.strandseq=conf[['most.frequent.state.strandseq']])
    
####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE
####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE
####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE ####### SWITCH TO PREVIOUS CODE FROM HERE
    
####### CHANGE FORMAT TO PREVIOUS ###########################################################################################
    
    # ID                -- same --
    # bins              -- same --
    # bincounts         has some extra columns ('state' and 'copy.number')  (artificially added)
    # qualityInfo       complexity is missing
    # segments          added 'num.bins' dropped 'GC' (This was an average of bins which does not make sense to me.)
    # weights           -- same --
    # distributions     -- same --
    # warnings          -- same --
    # breakpoints       to be added in next steps!!!
    
    model$bincounts     <- GRangesList()
    model$bincounts[[1]] <- model$bins
    
####### CONTINUE ############################################################################################################
    
    # Breakpoints and confidence intervals
    if (is.null(conf[['confint']])) {
      reads.file        <- NULL
    }else{
      reads.file        <- file.path(readspath, paste0(model$ID,'.RData'))
    }
    model$breakpoints <- getBreakpoints(model, fragments=reads.file, confint = conf[['confint']])
    model$qualityInfo <- as.list(getQC(model))
    save(model,file=file.path(path.method,basename(file.cur)))
  }



    modeldir <- file.path(modelpath, paste0('method-', method))
    if (conf[['refine.breakpoints']]) {
        refinedmodeldir <- file.path(refinedmodelpath, paste0('method-', method))
        if (!file.exists(refinedmodeldir)) { dir.create(refinedmodeldir, recursive=TRUE) }
    } else {
        refinedmodeldir <- modeldir
    }

    plotdir <- file.path(plotpath, paste0('method-', method))
    if (!file.exists(plotdir)) { dir.create(plotdir, recursive=TRUE) }
    
    browserdir <- file.path(browserpath, paste0('method-', method))
    if (!file.exists(browserdir)) { dir.create(browserdir, recursive=TRUE) }

    #========================
    ### refineBreakpoints ###
    #========================
    if (conf[['refine.breakpoints']]) {
        files <- list.files(modeldir, full.names=TRUE, pattern='.RData$')
        parallel.helper <- function(file) {
            tC <- tryCatch({
                savename <- file.path(refinedmodeldir,basename(file))
                if (!file.exists(savename)) {
                    model <- loadFromFiles(file)[[1]]
                    reads.file <- file.path(readspath, paste0(model$ID,'.RData'))
                    ## Refining breakpoints
                    if (!is.null(model$breakpoints)) {
                        message("Breakpoint refinement for ", basename(file))
                        model <- refineBreakpoints(model, fragments=reads.file, breakpoints=model$breakpoints, confint = conf[['confint']])
                    }
                    ptm <- startTimedMessage("Saving breakpoints to file ",savename," ...")
                    save(model, file=savename)
                    stopTimedMessage(ptm)
                }
            }, error = function(err) {
                stop(file,'\n',err)
            })
        }
        if (numcpu > 1) {
            ptm <- startTimedMessage("Refining breakpoints ...")
            temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
                parallel.helper(file)
            }
            stopTimedMessage(ptm)
        } else {
            # temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
            for (file in files) {
                parallel.helper(file)
            }
        }
    }
  
    #=======================
    ### Finding hotspots ###
    #=======================
    hslist <- list()
    if (is.null(hotspot.pval)) {
        message("Skipping breakpoint hotspot detection because 'hotspot.pval=NULL'.")
    } else {
        parallel.helper <- function(pattern) {
            ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
            ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
            breakpoints <- list()
            total.read.count <- numeric()
            for (file in ifiles) {
                hmm <- suppressMessages( loadFromFiles(file)[[1]] )
                breakpoints[[basename(file)]] <- hmm$breakpoints
                total.read.count[basename(file)] <- hmm$qualityInfo$total.read.count                                             # total.read.count is missing!!!!!
            }
            if (is.null(conf[['hotspot.bandwidth']])) {
                bw <- sum(as.numeric(seqlengths(hmm$bins))) / mean(total.read.count)
            } else {
                bw <- conf[['hotspot.bandwidth']][pattern]
            }
            hslist <- suppressMessages( hotspotter(breakpoints, bw=bw, pval=conf[['hotspot.pval']], spacing.bp = bw) )
            return(hslist)
        }
        if (numcpu > 1) {
            ptm <- startTimedMessage("Finding breakpoint hotspots ...")
            hslist <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
                parallel.helper(pattern)
            }
            stopTimedMessage(ptm)
        } else {
            ptm <- startTimedMessage("Finding breakpoint hotspots ...")
            hslist <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
                parallel.helper(pattern)
            }
            stopTimedMessage(ptm)
        }
        names(hslist) <- patterns
    }
  
    #===============
    ### Plotting ###
    #===============
    if (!file.exists(plotdir)) { dir.create(plotdir) }
    files <- list.files(refinedmodeldir, full.names=TRUE, pattern='.RData$')
    if (strandseq) {
        strandseq.string <- '_StrandSeq'
    } else {
        strandseq.string <- ''
    }
    if (conf[['refine.breakpoints']]) {
        strandseq.string <- paste0(strandseq.string, '_refined')
    }
  
    #------------------
    ## Plot heatmaps ##
    #------------------
    parallel.helper <- function(pattern) {
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        if (length(ifiles)>0) {
            savename=file.path(plotdir,paste0('genomeHeatmap_',sub('_$','',pattern), strandseq.string, '.pdf'))
            if (!file.exists(savename)) {
                suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.breakpoints=FALSE, hotspots=NULL, cluster=conf[['cluster.plots']]))
            }
        } else {
            warning("Plotting genomewide heatmaps: No files for pattern ",pattern," found.")
        }
    }
    if (numcpu > 1) {
        ptm <- startTimedMessage("Plotting genomewide heatmaps ...")
        temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
            parallel.helper(pattern)
        }
        stopTimedMessage(ptm)
    } else {
        # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
        for (pattern in patterns) {
            parallel.helper(pattern)
        }
    }
  
    parallel.helper <- function(pattern) {
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        if (length(ifiles)>0) {
            savename=file.path(plotdir,paste0('aneuploidyHeatmap_',sub('_$','',pattern), strandseq.string,'.pdf'))
            if (!file.exists(savename)) {
                grDevices::pdf(savename, width=30, height=max(0.3*length(ifiles), 2/2.54))
                ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=conf[['cluster.plots']]))
                print(ggplt)
                d <- grDevices::dev.off()
            }
        } else {
          warning("Plotting chromosome heatmaps: No files for pattern ",pattern," found.")
        }
    }
    if (numcpu > 1) {
        ptm <- startTimedMessage("Plotting chromosome heatmaps ...")
        temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
            parallel.helper(pattern)
        }
        stopTimedMessage(ptm)
    } else {
        ptm <- startTimedMessage("Plotting chromosome heatmaps ...")
        # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
        for (pattern in patterns) {
            parallel.helper(pattern)
        }
        stopTimedMessage(ptm)
    }
  
    #-------------------------
    ## Export browser files ##
    #-------------------------
    if (!file.exists(browserdir)) { dir.create(browserdir) }
    parallel.helper <- function(pattern) {
        ## Export CNV and breakpoints
        savename <- paste0(file.path(browserdir,sub('_$','',pattern)), strandseq.string)
        if (!file.exists(paste0(savename,'_CNV.bed.gz'))) {
            ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
            ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
            exportCNVs(ifiles, filename=savename, cluster=conf[['cluster.plots']], export.CNV=TRUE, export.breakpoints=TRUE)
        }
        ## Breakpoint hotspots
        savename <- file.path(browserdir,paste0(sub('_$','',pattern), strandseq.string, '_breakpoint-hotspots'))
        if (!file.exists(paste0(savename,'.bed.gz'))) {
            hotspots <- hslist[[pattern]]$hotspots
            if (!is.null(hotspots)) {
                exportGRanges(hotspots, filename=savename, trackname=basename(savename), score=hotspots$num.events, thickStart = hotspots$start.max, thickEnd = hotspots$end.max, priority=41)
            }
        }
        ## Hotspot densities
        savename <- file.path(browserdir,paste0(sub('_$','',pattern), strandseq.string, '_breakpoint-hotspot-densities'))
        if (!file.exists(paste0(savename,'.wig.gz'))) {
            densities <- hslist[[pattern]]$densities
            if (!is.null(hotspots)) {
                exportGRanges(densities, filename=savename, trackname=basename(savename), as.wiggle = TRUE, wiggle.val = densities$kde, priority=40)
            }
        }
    }
    if (numcpu > 1) {
      ptm <- startTimedMessage("Exporting browser files ...")
      temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
        parallel.helper(pattern)
      }
     stopTimedMessage(ptm)
    } else {
      # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
      for (pattern in patterns) {
        parallel.helper(pattern)
      }
    }
  
    #------------------
    ## Plot profiles ##
    #------------------
    parallel.helper <- function(pattern) {
      savename <- file.path(plotdir,paste0('profiles_',sub('_$','',pattern), strandseq.string,'.pdf'))
      if (!file.exists(savename)) {
        grDevices::pdf(file=savename, width=20, height=10)
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        for (ifile in ifiles) {
          tC <- tryCatch({
            model <- get(load(ifile))
            p1 <- graphics::plot(model, type='profile', plot.breakpoints=FALSE)
            p2 <- graphics::plot(model, type='histogram')
            cowplt <- cowplot::plot_grid(p1, p2, nrow=2, rel_heights=c(1.2,1))
            print(cowplt)
          }, error = function(err) {
            stop(ifile,'\n',err)
          })
        }
        d <- grDevices::dev.off()
      }
    }
    if (numcpu > 1) {
      ptm <- startTimedMessage("Making profile and distribution plots ...")
      temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
        parallel.helper(pattern)
      }
      stopTimedMessage(ptm)
    } else {
      ptm <- startTimedMessage("Making profile and distribution plots ...")
      # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
      for (pattern in patterns) {
        parallel.helper(pattern)
      }
      stopTimedMessage(ptm)
    }
  
    #--------------------
    ## Plot karyograms ##
    #--------------------
    parallel.helper <- function(pattern) {
      savename <- file.path(plotdir,paste0('karyograms_',sub('_$','',pattern), strandseq.string,'.pdf'))
      if (!file.exists(savename)) {
        grDevices::pdf(file=savename, width=12*1.4, height=2*4.6)
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        for (ifile in ifiles) {
          tC <- tryCatch({
            model <- get(load(ifile))
            print(graphics::plot(model, type='karyogram', plot.breakpoints=TRUE))
          }, error = function(err) {
            stop(ifile,'\n',err)
          })
        }
        d <- grDevices::dev.off()
      }
    }
    if (numcpu > 1) {
      ptm <- startTimedMessage("Plotting karyograms ...")
      temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
        parallel.helper(pattern)
      }
      stopTimedMessage(ptm)
    } else {
      ptm <- startTimedMessage("Plotting karyograms ...")
      # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
      for (pattern in patterns) {
        parallel.helper(pattern)
      }
      stopTimedMessage(ptm)
    }

  }


# ---------------------------------------------------------------------------------------------------------------------
# PART 3.02 | REFINE BREAKPOINTS
# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------
# PART 3.03 | FIND BREAKPOINT HOTSPOTS
# ---------------------------------------------------------------------------------------------------------------------



# =====================================================================================================================
# PART 4 OF 4 | CREATING BROWSER FILES, PLOTTING
# =====================================================================================================================

}                                                                                                                      # N: End of function.



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------