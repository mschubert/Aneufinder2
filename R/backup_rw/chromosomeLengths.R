chromosomeLengths <- function(assembly,chrom.type,input.files) {
    # N: You only indicate an assembly if you cannot get this information in
    # another way (from bam file). Also: A bam file does not need to have a
    # header section with chromosom length information!
    if (!is.null(assembly)) {
        if (is.character(assembly)) {
            # N: Chromosome length information from: file
            if (file.exists(assembly)) {
                dataf <- utils::read.table(assembly, sep='\t', header=TRUE)
                if(names(dataf)[1] != "chromosome" | names(dataf)[2] != "length")
                    stop("Incorrect column names assembly file or data frame. ",
                         "The names of the columns should be 'chromosome' and 'length'.")
            }
        } else {
            df.chroms <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(assembly)
            if (chrom.type == "chr")
                dataf <- df.chroms[,c('UCSC_seqlevel','UCSC_seqlength')]
            else
                dataf <- df.chroms[,c('NCBI_seqlevel','UCSC_seqlength')]
        }

    # N: Chromosome length information from: data frame
    } else if(is.data.frame(assembly)) {
        dataf <- assembly
        if(names(dataf)[1] != "chromosome" | names(dataf)[2] != "length")
            stop("Incorrect column names assembly file or data frame. The names ",
                 "of the columns should be 'chromosome' and 'length'.")
        else
            stop("'assembly' must be either a data.frame with columns ",
                 "'chromosome' and 'length' or a character specifying the ",
                 "assembly.")

        chrom.lengths <- dataf[,2]
        names(chrom.lengths) <- dataf[,1]

    # N: Chromosome length information from: bam file
    } else {
        bamfile <- grep('bam$', input.files, value=TRUE)[1]
        if(!is.na(bamfile))
            chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
        else
            stop("Not able to get chromosome length information.")
    }

    chrom.lengths
}
