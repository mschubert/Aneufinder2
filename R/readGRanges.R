#' Read a sequence region file into a GRanges object
#'
#' @param ...
readGRanges <- function(x, ..., seqinfo=GenomeInfoDb::seqinfo(x)) {
    UseMethod("readGRanges")
}

readGRanges.character <- function(x, ...) {
    if (length(x) != 1)
        stop("readGRanges needs exactly one file")

    if (grepl("\\.bam$", x))
        readGRanges(Rsamtools::BamFile(x))
    else if (grepl("\\.bed(\\.gz)?$", x))
        readGRanges.BedFile(x, ...)
    else
        stop("Unknown file extension (bam/bed supported): ", sQuote(x))
}

readGRanges.BamFile <- function(x, ..., seqinfo=GenomeInfoDb::seqinfo(x)) {
    if (!file.exists(paste0(x, ".bai"))) {
        ptm <- startTimedMessage("Couldn't find BAM index file, creating one.")
        Rsamtools::indexBam(x)
        stopTimedMessage(ptm)
    }

    ptm <- startTimedMessage("Reading file ", basename(x), " ...")
    args <- list(file=x, what="mapq", mapqFilter=min.mapq
                 param=Rsamtools::ScanBamParam(which=range(gr)))
    if (pairedEndReads)
        fun = GenomicAlignments::readGAlignmentPairs
    else
        fun = GenomicAlignments::readGAlignments
    if (remove.duplicate.reads)
        args$flag <- Rsamtools::scanBamFlag(isDuplicate=FALSE)
    data.raw = do.call(fun, args)
    stopTimedMessage(ptm)

    if (length(data.raw) == 0) {
        if (pairedEndReads)
            stop("No reads imported. Does your file really contain paired end ",
                 "reads? Try with 'pairedEndReads=FALSE'")
        else
            stop("No reads imported! Check your BAM-file ", x)
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

readGRanges.BedFile <- function(x, ..., seqinfo) {
    ptm <- startTimedMessage("Reading file ", basename(x), " ...")
    ccs <- c('character', 'numeric', 'numeric', 'NULL', 'integer', 'character')
    data.raw <- utils::read.table(bedfile, colClasses=ccs)
    data <- GenomicRanges::GRanges(
        seqnames = data.raw[,1],
        ranges = IRanges(start=data.raw[,2]+1, end=data.raw[,3]),
        strand = data.raw[,5],
        mcols = DataFrame(mapq = data.raw[,4]),
        seqinfo = seqinfo)
    stopTimedMessage(ptm)
    data
}

readGRanges.default <- function(x, ...) {
    stop("Could not read sequence information for type ", sQuote(class(x)))
}
