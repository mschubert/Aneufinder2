context("genome seqinfo")

test_that("genome information from bam file", {
    seqinfo = genome(bam)

    expect_setequal(names(seqinfo), c(1:19,'X'))
    expect_equal(GenomeInfoDb::seqlengths(seqinfo)[[1]], 195471971)
    expect_true(is.na(unique(GenomeInfoDb::genome(seqinfo)))) # not named
})

test_that("mm10", {
    seqinfo = genome("mm10")

    expect_equal(names(seqinfo), paste0("chr", c(1:19,'X')))
    expect_equal(GenomeInfoDb::seqlengths(seqinfo)[[1]], 195471971)
    expect_false(unique(GenomeInfoDb::isCircular(seqinfo)))
    expect_equal(unique(GenomeInfoDb::genome(seqinfo)), "mm10")
})
