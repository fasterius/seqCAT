library("CellAuthentication")
context("Read extracted variants")

# Read variants
file_1 <- system.file("extdata",
                      "extract.sample1.txt",
                      package = "CellAuthentication")
file_2 <- system.file("extdata",
                      "extract.sample2.txt",
                      package = "CellAuthentication")

granges_1 <- read_variants(file        = file_1,
                           sample_name = "sample1")

granges_2 <- read_variants(file        = file_2,
                           sample_name = "sample2")

# Tests
test_that("a GRanges object is returned", {
    expect_identical(class(granges_1)[1], "GRanges")
})

test_that("correct number of variants are read and de-duplicated", {
    expect_equal(length(granges_1), 383)
    expect_equal(length(granges_2), 382)
    expect_equal(length(mcols(granges_1)), 17)
    expect_equal(length(mcols(granges_2)), 17)
})

test_that("non-standard chromosomes are removed correcty", {
    expect_identical(levels(seqnames(granges_1)), "1")
    expect_identical(levels(seqnames(granges_2)), "1")
})
