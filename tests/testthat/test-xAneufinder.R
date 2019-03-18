context("Aneufinder main function")

test_that("integrated run", {
    skip("interactive test")

library(GenomicRanges)
library(GenomInfoDb)
source("helper-util.R")
source("../../R/Aneufinder.R")
source("../../R/rwConfig.R")
source("../../R/checkClass.R")
source("../../R/readGRanges.R")
source("../../R/genome.R")
source("../../R/partitionGenome.R")
source("../../R/binReads.R")
source("../../R/findCNVs.R")
source("../../R/findCNVs_edivisive.R")
source("../../R/util.R")
source("../../R/backup/collapseBins.R")
source("../../R/assign.distributions.R")
source("../../R/dnbinom.R")

    models = Aneufinder(bam, outputfolder=tempdir())
})
