context("parition genome")

test_that("fixed width bins", {
    genome = genome("mm10")
    bins = partitionGenome(genome, bin.size=1e6)

    expect_equal(length(bins), 2625)
    expect_equal(GenomeInfoDb::seqinfo(bins), genome)
    expect_setequal(as.character(GenomeInfoDb::seqnames(bins)), names(genome))
})

#test_that("variable width bins", {
#})

#test_that("error on genome mismatch", {
#})
