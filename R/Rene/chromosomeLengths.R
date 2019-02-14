
# ---------------------------------------------------------------------------------------------------------------------
#                                                  chromosomeLengths
# ---------------------------------------------------------------------------------------------------------------------
# Organization:  ERIBA (CRIPSR/iPSC facility)
# Programmer:    René Wardenaar (Original code written by Aaron Taudt)
# Starting date: 02-11-18
# Last modified: 02-11-18
# Version:       1.0
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ---------------------------------------------------------------------------------------------------------------------

RW_chromosomeLengths <- function(assembly,chrom.type,input.files){
  if(!is.null(conf[['assembly']])){                                                                                    # RW: You only indicate an assembly if you cannot get this information in another way (from bam file). Also: A bam file does not need to have a header section with chromosom length information!  
    if(is.character(conf[['assembly']])){
      if(file.exists(conf[['assembly']])){                                                                             # Chromosome length information from: file
        dataf <- utils::read.table(conf[['assembly']], sep='\t', header=TRUE)
        if(names(dataf)[1] != "chromosome" | names(dataf)[2] != "length"){
          stop("Incorrect column names assembly file or data frame. The names of the columns should be 'chromosome' and 'length'.")
        }
      }else{
        df.chroms <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(conf[['assembly']])
        if(chrom.type == "chr") {
          dataf <- df.chroms[,c('UCSC_seqlevel','UCSC_seqlength')]
        }else{
          dataf <- df.chroms[,c('NCBI_seqlevel','UCSC_seqlength')]
        }
      }
    }else if(is.data.frame(conf[['assembly']])){                                                                       # Chromosome length information from: data frame
      dataf <- conf[['assembly']]
      if(names(dataf)[1] != "chromosome" | names(dataf)[2] != "length"){
        stop("Incorrect column names assembly file or data frame. The names of the columns should be 'chromosome' and 'length'.")
      }
    }else{
      stop("'assembly' must be either a data.frame with columns 'chromosome' and 'length' or a character specifying the assembly.")
    }
    chrom.lengths        <- dataf[,2]
    names(chrom.lengths) <- dataf[,1]
  }else{
    bamfile <- grep('bam$', input.files, value=TRUE)[1]
    if(!is.na(bamfile)){                                                                                               # Chromosome length information from: bam file
      chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
    }else{
      stop("Not able to get chromosome length information.")
    }
  }
  return(chrom.lengths)
}



# ---------------------------------------------------------------------------------------------------------------------
# END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
# ---------------------------------------------------------------------------------------------------------------------