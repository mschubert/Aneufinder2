fixedWidthBins <- function(bamfile=NULL, assembly=NULL, chrom.lengths=NULL, chromosome.format, binsizes=1e6, stepsizes=NULL, chromosomes=NULL) {

	### Check user input ###
  if (length(binsizes) == 0) {
    return(list())
  }
	if (is.null(bamfile) & is.null(assembly) & is.null(chrom.lengths)) {
		stop("Please specify either a 'bamfile', 'assembly' or 'chrom.lengths'")
	}
	if (is.null(bamfile) & is.null(chrom.lengths)) {
		trigger.error <- chromosome.format
	}
  if (!is.null(stepsizes)) {
    if (length(stepsizes) != length(binsizes)) {
      stop("Need one element in 'stepsizes' for each element in 'binsizes'.")
    }
    if (any(binsizes < stepsizes)) {
        stop("'stepsizes' must be smaller/equal than 'binsizes'")
    }
  }

	### Get chromosome lengths ###
	if (!is.null(bamfile)) {
      ptm <- startTimedMessage(paste0("Reading header from ", bamfile, " ..."))
      chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
      stopTimedMessage(ptm)
	} else if (!is.null(assembly)) {
		if (is.character(assembly)) {
			ptm <- startTimedMessage("Fetching chromosome lengths from UCSC ...")
			df <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(assembly)
			stopTimedMessage(ptm)
		} else if (is.data.frame(assembly)) {
			df <- assembly
		} else {
			stop("Unknown assembly")
		}
		chrom.lengths <- df$UCSC_seqlength
		if (chromosome.format=='UCSC') {
			names(chrom.lengths) <- df$UCSC_seqlevel
		} else if (chromosome.format=='NCBI') {
			names(chrom.lengths) <- df$NCBI_seqlevel
		}
		chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
		chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
	} else if (!is.null(chrom.lengths)) {
		chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
		chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
	}
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if none of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('Could not find length information for any of the specified chromosomes: ', chrstring, '. Pay attention to the naming convention in your data, e.g. "chr1" or "1".')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning('Could not find length information for the following chromosomes: ', diffs)
	}

	### Making fixed-width bins ###
	bins.list <- list()
	for (ibinsize in 1:length(binsizes)) {
	  binsize <- binsizes[ibinsize]
    ptm <- startTimedMessage("Making fixed-width bins for bin size ", binsize, " ...")
    chrom.lengths.floor <- floor(chrom.lengths / binsize) * binsize
    bins <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor[chroms2use], tilewidth=binsize), use.names=FALSE)
    bins <- bins[end(bins) > 0] # no chromosomes that are smaller than binsize
    if (any(width(bins)!=binsize)) {
      stop("tileGenome failed")
    }
		# seqlengths(bins) <- as.integer(chrom.lengths[names(seqlengths(bins))])
    seqlengths(bins) <- chrom.lengths[chroms2use]
    if (!is.null(stepsizes)) {
      shift.bp <- 0
  	  stepsize <- stepsizes[ibinsize]
  	  bins.list.step <- GRangesList()
  	  while (shift.bp < binsize) {
  	    bins.list.step[[as.character(shift.bp)]] <- suppressWarnings( trim(shift(bins, shift.bp)) )
        shift.bp <- stepsize + shift.bp
  	  }
  	  bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE), '_stepsize_', format(stepsize, scientific=TRUE, trim=TRUE))]] <- bins.list.step
    } else {
  		bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE))]] <- bins
    }

    skipped.chroms <- setdiff(seqlevels(bins), as.character(unique(seqnames(bins))))
    if (length(skipped.chroms)>0) {
        warning("The following chromosomes were skipped because they are smaller than binsize ", binsize, ": ", paste0(skipped.chroms, collapse=', '))
    }
    stopTimedMessage(ptm)

	}

	return(bins.list)

}


variableWidthBins <- function(reads, binsizes, stepsizes=NULL, chromosomes=NULL) {
	
	### Check user input ###
	chroms.in.data <- seqlevels(reads)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if non of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('Could not find length information for any of the specified chromosomes: ', chrstring)
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning('Could not find length information for the following chromosomes: ', diffs)
	}
  if (!is.null(stepsizes)) {
    if (length(stepsizes) != length(binsizes)) {
      stop("Need one element in 'stepsizes' for each element in 'binsizes'.")
    }
    if (any(binsizes < stepsizes)) {
        stop("'stepsizes' must be smaller/equal than 'binsizes'")
    }
  }

	## Drop unwanted seqlevels
	reads <- reads[seqnames(reads) %in% chroms2use]
	reads <- keepSeqlevels(reads, chroms2use)

	## Make fixed width bins
	ptm <- startTimedMessage("Binning reads in fixed-width windows ...")
	binned.list <- suppressMessages( binReads(reads, assembly=NULL, binsizes=binsizes, stepsizes=stepsizes, calc.complexity=FALSE, chromosomes=chromosomes) )
	stopTimedMessage(ptm)
	
	## Sort the reads
	strand(reads) <- '*'
	reads <- sort(reads)
	
	## Loop over binsizes
	bins.list <- list()
	for (i1 in 1:length(binsizes)) {
		binsize <- binsizes[i1]
		if (is.null(stepsizes)) {
  		ptm <- startTimedMessage("Making variable-width windows for bin size ", binsize, " ...")
		} else {
  		ptm <- startTimedMessage("Making variable-width windows for bin size ", binsize, " and step size ", stepsizes[i1], " ...")
		}
		if (is(binned.list[[i1]], "GRangesList")) {
  		binned <- binned.list[[i1]][[1]]
		} else if (is(binned.list[[i1]], "GRanges")) {
  		binned <- binned.list[[i1]]
		}
		## Get median of histogram
		mediancount <- as.integer(median(binned$counts[binned$counts>0]))
		mediancount.perstep <- 0
		if (!is.null(stepsizes)) {
		  stepsize <- stepsizes[i1]
		  numsteps <- as.integer(binsize / stepsize)
  		mediancount.perstep <- unique(as.integer((0:(numsteps-1)) * mediancount / numsteps))
		}
		bins.list.step <- GRangesList()
		for (istep in 1:length(mediancount.perstep)) {
  		## Pick only every mediancount read
  		subreads <- GRangesList()
  		skipped.chroms <- character()
  		for (chrom in chroms2use) {
  			reads.chr <- reads[seqnames(reads)==chrom]
  			if (length(reads.chr) >= mediancount) {
  			  if (mediancount.perstep[istep] == 0) {
    				idx <- seq(mediancount, length(reads.chr), by=mediancount)
  			  } else {
    				idx <- seq(mediancount.perstep[istep], length(reads.chr), by=mediancount)
  			  }
  				subreads[[chrom]] <- reads.chr[idx]
  			} else {
  				skipped.chroms[chrom] <- chrom
  			}
  		}
  		subreads <- unlist(subreads, use.names=FALSE)
  		## Adjust length of reads to get consecutive bins
  		subreads <- resize(subreads, width=1)
  		## Make new bins
  		bins <- gaps(subreads, start=1L, end=seqlengths(subreads)-1L) # gaps until seqlengths-1 because we have to add 1 later to get consecutive bins
  		bins <- bins[strand(bins)=='*']
  		end(bins) <- end(bins) + 1
  		## Remove first bin
  		if (mediancount.perstep[istep] > 0) {
  		  bins <- bins[-1]
  		}
  		## We don't want incomplete bins at the end
  		bins.split <- split(bins, seqnames(bins))
  		bins.split <- endoapply(bins.split, function(x) { x[-length(x)] })
  		bins <- unlist(bins.split, use.names=FALSE)
  		## Remove skipped chromosomes
  		bins <- bins[!seqnames(bins) %in% skipped.chroms]
  		bins <- keepSeqlevels(bins, setdiff(seqlevels(bins), skipped.chroms))
  		bins.list.step[[as.character(istep)]] <- bins
		}
		if (length(skipped.chroms)>0) {
			warning("The following chromosomes were skipped because they are smaller than binsize ", binsize, ": ", paste0(skipped.chroms, collapse=', '))
		}

		if (is.null(stepsizes)) {
  		bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE))]] <- bins.list.step[[1]]
		} else {
  		bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE), '_stepsize_', format(stepsize, scientific=TRUE, trim=TRUE))]] <- bins.list.step
		}
		stopTimedMessage(ptm)
	}
	
	return(bins.list)

}


