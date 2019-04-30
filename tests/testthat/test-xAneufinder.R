context("Aneufinder main function")

test_that("integrated run", {
    skip("interactive test")

#library(AneuFinder2)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
source("helper-util.R")
source("../../R/Aneufinder.R")
source("../../R/rwConfig.R")
source("../../R/checkClass.R")
source("../../R/readGRanges.R")
source("../../R/genome.R")
source("../../R/partitionGenome.R")
source("../../R/addGCcontent.R")
source("../../R/correctGC.R")
source("../../R/binReads.R")
source("../../R/findCNVs.R")
source("../../R/findCNVs_edivisive.R")
source("../../R/findCNVs_HMM.R")
source("../../R/util.R")
source("../../R/distributions.R")
source("../../R/util.R")
source("../../R/genomeBins.R")
source("../../R/adjustBins.R")
source("../../R/collapseBins.R") #TODO: rewrite
#source("../../R/backup/clusterHMMs.R")
#source("../../R/backup/qualityMeasures.R")
source("../../R/plotKaryogram.R") #TODO: remove this + replace by commented code in plotKaryograms
source("../../R/plotKaryograms.R")
source("../../R/plotProfile.R")

    dir = tempdir()
    models = Aneufinder(bam, assembly="GRCm38", outputfolder=dir,
                        variable.width.reference=bam)
})
