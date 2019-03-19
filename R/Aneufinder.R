Aneufinder <- function(inputfolder, outputfolder, configfile=NULL, numCPU=1,
                       reuse.existing.files=TRUE, binsizes=1e6,
                       stepsizes=binsizes, variable.width.reference=NULL,
                       reads.per.bin=NULL, pairedEndReads=FALSE, assembly=NULL,
                       chromosomes=NULL, remove.duplicate.reads=TRUE,
                       min.mapq=10, blacklist=NULL, use.bamsignals=FALSE,
                       reads.store=FALSE, correction.method=NULL,
                       GC.BSgenome=NULL, method='edivisive', strandseq=FALSE,
                       R=10, sig.lvl=0.1, eps=0.01, max.time=60, max.iter=5000,
                       num.trials=15,
                       states=c('zero-inflation',paste0(0:10,'-somy')),
                       most.frequent.state='2-somy',
                       most.frequent.state.strandseq='1-somy', confint=NULL,
                       refine.breakpoints=FALSE, hotspot.bandwidth=NULL,
                       hotspot.pval=5e-2, cluster.plots=TRUE) {

    # load config params, but arguments take priority
    conf <- readConfig(configfile)
    args <- lapply(as.list(match.call())[-1], eval, envir=parent.frame())
    conf <- utils::modifyList(conf, args)
    conf$reuse.existing.files <- reuse.existing.files #FIXME: take this from call
    conf$cluster.plots <- cluster.plots
#    checkClass(conf=conf) # this is still too verbose
    # ^^ also should check/set [most will be done in checkClass already]:
    #  binsize/stepsize is integer >= 1; one step size for each bin size
    #  reads.per.bin >= 1, integer
    #  min.mapq >= 0, integer
    #  correction.method is NULL or "GC"; BSgenome is available (fail early)
    #  reuse.existing.directory should only delete aneufinder files, stop if others
    #  method %in% HMM, dnacopy, edivisive
    #  R < 1 (edivisive param)
    #  0 < sig.lvl <= 1
    #  0 < eps < 0.1
    #  max.time -1 or >=1, is integer
    #  num.trials integer > 0
    #  states must have zero-inflation [first], [0-9]+-somy
    #  set state.num (numeric states) that must be ordered [not incl state'+']
    #  most.frequent.state[.strandseq] must be in states
    #  warn if chromosome requested but not in seq file (error if no matches)
    #  NCBI/UCSC chromosome names consistent
    #  blacklist also has right chromosomes

    makedir(conf$outputfolder)

    if (is.null(assembly)) # this currently doesn't work with a inputFOLDER and no assembly
        seqinfo <- genome(inputfolder[1])
    else {
        seqinfo <- genome(assembly)
        assembly <- unique(GenomeInfoDb::genome(seqinfo))
    }

    ###
    ### Bin the genome
    ###
    fname <- args2fname(file.path(conf$outputfolder, assembly),
        binsize=binsizes, stepsize=stepsizes)
    if (file.exists(fname) && conf$reuse.existing.files) {
        bins <- readRDS(fname)
    } else {
        bins <- partitionGenome(seqinfo, binsize=binsizes,
                                reads.per.bin=reads.per.bin, stepsize=stepsizes)
        if ("GC" %in% correction.method)
            bins <- addGCcontent(bins, BSgenome=GC.BSgenome)
        saveRDS(bins, file=fname)
    }

    # missing: filtered reads -- file.path(conf[['outputfolder']],'filtered')

    ###
    ### Assign reads to bins
    ###
    bin_reads <- function(x, dir=makedir(conf$outputfolder, "binned")) {
        fname <- args2fname(file.path(dir, basename(x)),
            binsize=binsizes, stepsize=stepsizes)
        if (file.exists(fname) && conf$reuse.existing.files) {
            reads <- readRDS(fname)
        } else {
            args <- conf[intersect(names(conf), names(formals(readGRanges)))]
            reads <- do.call("binReads", c(args, list(reads=x, bins=bins)))
            if ("GC" %in% correction.method)
                reads <- correctGC(reads, method="loess")
            saveRDS(reads, file=fname)
        }
    }
    reads <- lapply(inputfolder, bin_reads)

    ###
    ### Identify CNVs
    ###
    find_cnvs <- function(x, dir=makedir(conf$outputfolder, "MODELS")) {
        fname <- args2fname(file.path(dir, attr(x, 'ID')), method=method)
        if (file.exists(fname) && conf$reuse.existing.files) {
            models <- readRDS(fname)
        } else {
            args <- conf[intersect(names(conf), names(formals(findCNVs)))]
            models <- do.call("findCNVs", c(args, list(binned=x)))
            saveRDS(models, file=fname)
        }
    }
    models <- lapply(reads, find_cnvs) #FIXME: fails on first run

    ###
    ### Refine breakpoints, breakpoint hotspots
    ###

    ###
    ### Make plots
    ###
    plotdir <- makedir(conf$outputfolder, "PLOTS")
    fname <- args2fname(file.path(plotdir, "genomeHeatmap"), ext=".pdf")
    if (!file.exists(fname))
        heatmapGenomewide(models, cluster=conf$cluster.plots, file=fname)

    fname <- args2fname(file.path(plotdir, "aneuploidyHeatmap"), ext=".pdf")
    if (!file.exists(fname)) {
        pdf(fname, width=30, height=max(0.3*length(models), 2/2.54))
        for (i in seq_along(models)) {
            message("Read density plot for: ", names(models)[i])
            print(plot(models[[i]]))
        }
        dev.off()
    }
}
