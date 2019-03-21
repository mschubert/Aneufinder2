#' Read a sequence region file into a GRanges object
#'
#' @param x  File name to read into GRanges object
#' @param seqinfo  Seqinfo object with metadata about the genome
#' @return  A GRanges object containing all the reads of the sequence file
readGRanges <- function(x, seqinfo=GenomeInfoDb::seqinfo(x), ...) {
    UseMethod("readGRanges")
}

readGRanges.character <- function(x, ...) {
    if (length(x) != 1)
        stop("readGRanges needs exactly one file")

    if (grepl("\\.bam$", x))
        gr <- readGRanges(Rsamtools::BamFile(x))
    else if (grepl("\\.bed(\\.gz)?$", x))
        gr <- readGRanges.BedFile(x, ...)
    else
        stop("Unknown file extension (bam/bed supported): ", sQuote(x))

    attr(gr, 'ID') <- tools::file_path_sans_ext(basename(x))
    gr
}

readGRanges.BamFile <- function(x, seqinfo=GenomeInfoDb::seqinfo(x), min.mapq=10,
                                pairedEndReads=FALSE, remove.duplicate.reads=FALSE,
                                max.fragment.width=1000) {
    if (!file.exists(paste0(x$path, ".bai"))) { # x$index?
        ptm <- startTimedMessage("Couldn't find BAM index file, creating one.")
        Rsamtools::indexBam(x)
        stopTimedMessage(ptm)
    }

    ptm <- startTimedMessage("Reading file ", basename(x$path), " ...")
    args <- list(file=x) #, what="mapq", mapqFilter=min.mapq),
                 #param=Rsamtools::ScanBamParam(which=range(gr)))
    if (pairedEndReads)
        fun = GenomicAlignments::readGAlignmentPairs
    else
        fun = GenomicAlignments::readGAlignments
#    if (remove.duplicate.reads)
#        args$flag <- Rsamtools::scanBamFlag(isDuplicate=FALSE)
    data.raw = do.call(fun, args)
    stopTimedMessage(ptm)

    if (length(data.raw) == 0) {
        if (pairedEndReads)
            stop("No reads imported. Does your file really contain paired end ",
                 "reads? Try with 'pairedEndReads=FALSE'")
        else
            stop("No reads imported! Check your BAM-file ", basename(x$path))
    }

    ptm <- startTimedMessage("Converting to GRanges ...")
    if (pairedEndReads) # treat as one fragment
        data <- GenomicAlignments::granges(data.raw, use.mcols=TRUE,
                                           on.discordant.seqnames='drop')
    else
        data <- GenomicAlignments::granges(data.raw, use.mcols=TRUE)
    stopTimedMessage(ptm)

    ptm <- startTimedMessage("Filtering reads ...")
    data <- data[GenomicRanges::width(data) <= max.fragment.width]
    stopTimedMessage(ptm)

    data
}

readGRanges.BedFile <- function(x, seqinfo=NULL) {
    ptm <- startTimedMessage("Reading file ", basename(x), " ...")
    ccs <- c('character', 'numeric', 'numeric', 'NULL', 'integer', 'character')
    data.raw <- utils::read.table(x, colClasses=ccs)
    data <- GenomicRanges::GRanges(
        seqnames = data.raw[,1],
        ranges = IRanges::IRanges(start=data.raw[,2]+1, end=data.raw[,3]),
        strand = data.raw[,5],
        mcols = S4Vectors::DataFrame(mapq = data.raw[,4]),
        seqinfo = seqinfo)
    stopTimedMessage(ptm)
    data
}

readGRanges.default <- function(x, ...) {
    stop("Could not read sequence information for type ", sQuote(class(x)))
}
