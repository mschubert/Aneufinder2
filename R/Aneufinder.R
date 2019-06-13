#' Wrapper function for the \code{\link{AneuFinder}} package
#'
#' This function is an easy-to-use wrapper to \link[AneuFinder:binning]{bin the
#' data}, \link[AneuFinder:findCNVs]{find copy-number-variations},
#' \link[AneuFinder:getBreakpoints]{locate breakpoints}, plot
#' \link[AneuFinder:heatmapGenomewide]{genomewide heatmaps},
#' \link[AneuFinder:plot.aneuHMM]{distributions, profiles and karyograms}.
#'
#' @param inputfolder Folder with either BAM or BED files.
#' @param outputfolder Folder to output the results. If it does not exist it
#'   will be created.
#' @param configfile A file specifying the parameters of this function (without
#'   \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the
#'   parameters in a file can be handy if many samples with the same parameter
#'   settings are to be run. If a \code{configfile} is specified, it will take
#'   priority over the command line parameters.
#' @param numCPU The numbers of CPUs that are used. Should not be more than
#'   available on your machine.
#' @param reuse.existing.files A logical indicating whether or not existing
#'   files in \code{outputfolder} should be reused.
#' @inheritParams readGRanges
#' @inheritParams binReads
#' @param reads.store Set \code{reads.store=TRUE} to
#'   store read fragments as RData in folder 'data' and as BED files in
#'   'BROWSERFILES/data'. This option will force \code{use.bamsignals=FALSE}.
#' @param correction.method Correction methods to be used for the binned read
#'   counts. Currently only \code{'GC'}.
#' @param GC.BSgenome A \code{BSgenome}
#'   object which contains the DNA sequence that is used for the GC correction.
#' @param strandseq A logical indicating whether the data comes from Strand-seq
#'   experiments. If \code{TRUE}, both strands carry information and are treated
#'   separately.
#' @inheritParams edivisive.findCNVs
#' @inheritParams HMM.findCNVs
#' @inheritParams findCNVs
#' @param confint Desired confidence interval for breakpoints. Set
#'   \code{confint=NULL} to disable confidence interval estimation. Confidence
#'   interval estimation will force \code{reads.store=TRUE}.
#' @param refine.breakpoints A logical indicating whether breakpoints from the
#'   HMM should be refined with read-level information.
#'   \code{refine.breakpoints=TRUE} will force \code{reads.store=TRUE}.
#' @param hotspot.bandwidth A vector the same length as \code{binsizes} with
#'   bandwidths for breakpoint hotspot detection (see \code{\link{hotspotter}}
#'   for further details). If \code{NULL}, the bandwidth will be chosen
#'   automatically as the average distance between reads.
#' @param hotspot.pval
#'   P-value for breakpoint hotspot detection (see \code{\link{hotspotter}} for
#'   further details). Set \code{hotspot.pval = NULL} to skip hotspot detection.
#' @param cluster.plots A logical indicating whether plots should be clustered
#'   by similarity.
#' @return \code{NULL}
#' @author Aaron Taudt
#' @import foreach
#' @import doParallel
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot
#' @importFrom utils read.table write.table
#' @importFrom cowplot plot_grid
#' @export
#'
#' @examples
#' \dontrun{
#' ## The following call produces plots and genome browser files for all BAM files in "my-data-folder"
#' Aneufinder(inputfolder="my-data-folder", outputfolder="my-output-folder")
#' }
Aneufinder <- function(inputfolder, outputfolder, assembly, configfile=NULL,
                       numCPU=1, reuse.existing.files=TRUE, binsizes=1e6,
                       stepsizes=binsizes, variable.width.reference=NULL,
                       reads.per.bin=NULL, pairedEndReads=FALSE,
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
    conf <- list()
    if (!is.null(configfile))
        conf <- RcppTOML::parseTOML(configfile)
    args <- lapply(as.list(match.call())[-1], eval, envir=parent.frame())
    conf <- utils::modifyList(conf, args)
    conf$reuse.existing.files <- reuse.existing.files #FIXME: take this from call (match w/ default args)
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
    seqinfo <- genome(assembly)

    if (length(inputfolder) == 1 && dir.exists(inputfolder))
        inputfolder = list.files(inputfolder, "\\.(bam|bed(\\.gz)?)$", full.names=TRUE)
    names(inputfolder) <- tools::file_path_sans_ext(basename(inputfolder))

    ###
    ### Bin the genome
    ###
    fname <- args2fname(file.path(conf$outputfolder, assembly),
        "fixed", binsize=binsizes, stepsize=stepsizes)
    if (file.exists(fname) && conf$reuse.existing.files) {
        bins <- readRDS(fname)
        changed <- FALSE
    } else {
        bins <- genomeBins(seqinfo, binsize=binsizes,
                           reads.per.bin=reads.per.bin, stepsize=stepsizes)
        changed <- TRUE
    }
    if ("GC" %in% correction.method && is.null(variable.width.reference) &&
            ! "GC.content" %in% colnames(bins)) {
        bins <- addGCcontent(bins, BSgenome=GC.BSgenome)
        changed <- TRUE
    }
    if (changed)
        saveRDS(bins, file=fname)

    if (!is.null(variable.width.reference)) {
        reads <- readGRanges(variable.width.reference)
        fname <- args2fname(file.path(conf$outputfolder, assembly),
            "variable", binsize=binsizes, stepsize=stepsizes)
        if (file.exists(fname) && conf$reuse.existing.files) {
            bins <- readRDS(fname)
            changed <- FALSE
        } else {
            bins <- adjustBins(bins, reads)
            changed <- TRUE
        }
        if ("GC" %in% correction.method && ! "GC.content" %in% colnames(bins)) {
            bins <- addGCcontent(bins, BSgenome=GC.BSgenome)
            changed <- TRUE
        }
        if (changed)
            saveRDS(bins, file=fname)
    }

    # missing: filtered reads -- file.path(conf[['outputfolder']],'filtered')
    if (is.null(chromosomes))
        bins = GenomeInfoDb::keepStandardChromosomes(bins, pruning.mode="coarse")
    else
        bins = GenomeInfoDb::keepSeqlevels(bins, chromosomes, pruning.mode="coarse")

    ###
    ### Assign reads to bins
    ###
    bin_reads <- function(x, dir=makedir(conf$outputfolder, "binned")) {
        fname <- args2fname(file.path(dir, tools::file_path_sans_ext(basename(x))),
            binsize=binsizes, stepsize=stepsizes)
        if (file.exists(fname) && conf$reuse.existing.files)
            return(readRDS(fname))

        args <- conf[intersect(names(conf), names(formals(readGRanges)))]
        reads <- do.call("binReads", c(args, list(reads=x, bins=bins)))
        if ("GC" %in% correction.method)
            reads <- correctGC(reads, method="loess")
        saveRDS(reads, file=fname)
        reads
    }
    reads <- lapply(inputfolder, bin_reads)

    ###
    ### Identify CNVs
    ###
    find_cnvs <- function(x, dir=makedir(conf$outputfolder, "MODELS")) {
        fname <- args2fname(file.path(dir, attr(x, 'ID')), method=method)
        if (file.exists(fname) && conf$reuse.existing.files)
            return(readRDS(fname))

        args <- conf[intersect(names(conf), names(formals(findCNVs)))]
        models <- do.call("findCNVs", c(args, list(binned=x)))
        saveRDS(models, file=fname)
        models
    }
    models <- lapply(reads, find_cnvs)

    ###
    ### Refine breakpoints, breakpoint hotspots
    ###

    ###
    ### Make plots
    ###
    plotdir <- makedir(conf$outputfolder, "PLOTS")

    fname <- args2fname(file.path(plotdir, "karyograms"), ext=".pdf")
    if (!file.exists(fname)) {
        pdf(fname, width=20, height=length(models)+4)
        #print(plotKaryograms(models, cluster=conf$cluster.plots))
        for(i in seq_along(models)){
            print(plotKaryogram(models[[i]]))
        }
        dev.off()
    }

    fname <- args2fname(file.path(plotdir, "profile"), ext=".pdf")
    if (!file.exists(fname)) {
        pdf(fname, width=30, height=max(0.3*length(models), 2/2.54))
        for (i in seq_along(models)) {
            message("Read density plot for: ", names(models)[i])
            print(plotProfile(models[[i]]))
        }
        dev.off()
    }
}
