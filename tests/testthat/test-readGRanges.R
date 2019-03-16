context("reading files")

fdir = system.file("extdata", package="AneuFinderData")
bam = file.path(fdir, "BB150803_IV_074.bam")
bed = file.path(fdir, "KK150311_VI_07.bam.bed.gz")

test_that("reading bam file into GRanges", {
    ranges = readGRanges(bam)
    expect_equal(c(class(ranges)), "GRanges")
})

test_that("reading bed file into GRanges", {
    ranges = readGRanges(bed)
    expect_equal(c(class(ranges)), "GRanges")
})

test_that("error on genome mismatch", {
})
