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

    if (length(inputfolder) == 1 && dir.exists(inputfolder))
        inputfolder <- list.files(inputfolder, "\\.(bam|bed(\\.gz)?)$", full.names=TRUE)
    if (length(inputfolder) == 0)
        stop("No BAM or BED files supplied or in directory")
    names(inputfolder) <- tools::file_path_sans_ext(basename(inputfolder))

    if (!file.exists(conf[['outputfolder']]))
        dir.create(conf[['outputfolder']])

    if (is.null(assembly))
        seqinfo <- genome(inputfolder[1])
    else
        seqinfo <- genome(assembly)

    # provide a function to call(partitionGenome, conf) that matches args?
    #   also potentially combine with looking if results are there + skip (otherwise +save)
    bins <- partitionGenome(seqinfo, binsize=binsizes,
                            reads.per.bin=reads.per.bin, stepsize=stepsizes)
    if ("GC" %in% correction.method)
        bins <- addGCcontent(bins, BSgenome=GC.BSgenome)

    # check if the binned files are available first and load, or save otherwise
    # probably provide a filenames.S3 to query file names based on function calls
    #   file.path(conf[['outputfolder']],'filtered')
    #   save(reads,file=file.path(path.filtered.reads,paste0(basename(file.cur),'.Rdata')))
    #   path.uncorrected.bins   <- file.path(conf[['outputfolder']],'binned')
    #   paste0("binsize_",format(binsize,scientific=TRUE,trim=TRUE),"_stepsize_",
    #                            format(stepsize,scientific=TRUE,trim=TRUE))
    #  inp_file            <- basename(file.cur)
    #  inp_file            <- substr(inp_file,1,(nchar(inp_file)-6))
    #  file.save           <- file.path(path.uncorrected.bins,paste0(inp_file,"_",combi,".RData"))
    args <- conf[intersect(names(conf), names(formals(readGRanges)))]
    reads <- do.call("binReads", c(args, list(reads=inputfolder, bins=bins)))

    if ("GC" %in% correction.method)
        reads <- correctGC(reads, method="loess")

    # reads: {output}/data
    # create {output}/MODELS{,_refined} /PLOTS /BROWSERFILES
    # save

    args <- conf[intersect(names(conf), names(formals(findCNVs)))]
    models <- do.call("findCNVs", c(args, list(binned=reads)))

    # missing: refine breakpoints, breakpoint hotspots

    # create plotdir

    fname <- file.path(plotdir,paste0('genomeHeatmap_',sub('_$','',pattern), '.pdf'))
    heatmapGenomewide(models, cluster=TRUE, file=fname)

    fname <- file.path(plotdir,paste0('aneuploidyHeatmap_',sub('_$','',pattern), '.pdf'))
    pdf(fname, width=30, height=max(0.3*length(models), 2/2.54))
    for (i in seq_along(all_models)) {
        message("Read density plot for: ", names(all_models)[i])
        if (class(models[[i]]) == "aneuHMM") #TODO: is this required?
            print(plot(models[[i]]))
    }
    dev.off()
}
