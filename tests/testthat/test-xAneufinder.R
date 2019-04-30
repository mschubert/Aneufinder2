context("Aneufinder main function")

test_that("integrated run", {
    skip("interactive test")

library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
source("helper-util.R")
for (fname in list.files("../../R", "\\.R", full.names=TRUE))
    source(fname)
#source("../../R/collapseBins.R") #TODO: rewrite
#source("../../R/backup/clusterHMMs.R")
#source("../../R/backup/qualityMeasures.R")

    dir = tempdir()
    models = Aneufinder(bam, assembly="GRCm38", outputfolder=dir,
                        variable.width.reference=bam)
})
