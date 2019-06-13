This folder contains the following files:
-----------------------------------------
- chrominfo.tsv: A tab-separated file with chromosome lengths.
- AneuFinder.config: A text file with all the parameters that were used to run Aneufinder().

This folder contains the following folders. Some folders are [optional] depending on the parameters:
-------------------------------------------
- binned: RData files with the results of the binnig step. Contains GRanges objects with binned genomic coordinates and read counts.
- [binned-GC]: RData files with the results of the GC-correction step. Depends on option 'correction.method="GC"'. Contains GRanges objects with binned genomic coordinates and read counts.
- BROWSERFILES: Bed files for upload to the UCSC genome browser.
- [data]: RData files with GRanges containing read-fragments. Depends on option 'reads.store=TRUE'.
- MODELS [or MODELS-StrandSeq]: RData files with aneuHMM objects. Result of the copy-number and breakpoint estimation step. Depends on option 'strandseq=TRUE'.
- MODELS_refined [or MODELS-StrandSeq_refined]: RData files with aneuHMM objects after the breakpoint refinement step. Depends on option 'strandseq=TRUE'.
- PLOTS: Several plots that are produced by default.
