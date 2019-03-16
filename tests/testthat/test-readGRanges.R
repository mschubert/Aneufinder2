context("reading files")

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
