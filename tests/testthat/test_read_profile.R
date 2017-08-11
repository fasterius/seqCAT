library("CellAuthentication")
context("Read SNV profiles")

# Read variants
file_1 <- system.file("extdata",
                      "profile_1.txt.gz",
                      package = "CellAuthentication")
file_2 <- system.file("extdata",
                      "profile_2.txt.gz",
                      package = "CellAuthentication")

profile_1 <- read_profile(file        = file_1,
                          sample_name = "sample1")

profile_2 <- read_profile(file        = file_2,
                          sample_name = "sample2")

# Tests
test_that("a GRanges object is returned", {
    expect_identical(class(profile_1)[1], "GRanges")
})

test_that("correct number of variants are read and de-duplicated", {
    expect_equal(length(profile_1), 383)
    expect_equal(length(profile_2), 382)
    expect_equal(length(mcols(profile_1)), 17)
    expect_equal(length(mcols(profile_2)), 17)
})

test_that("non-standard chromosomes are removed correcty", {
    expect_identical(levels(seqnames(profile_1)), "1")
    expect_identical(levels(seqnames(profile_2)), "1")
})
