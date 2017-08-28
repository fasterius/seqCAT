library("CellAuthentication")
context("Read SNV profiles")

# Read variants
file_1 <- system.file("extdata",
                      "test_profile_1.txt.gz",
                      package = "CellAuthentication")
file_2 <- system.file("extdata",
                      "test_profile_2.txt.gz",
                      package = "CellAuthentication")

test_profile_1 <- read_profile(file        = file_1,
                               sample_name = "sample1")

test_profile_2 <- read_profile(file        = file_2,
                               sample_name = "sample2")

# Tests
test_that("a GRanges object is returned", {
    expect_identical(class(test_profile_1)[1], "GRanges")
})

test_that("correct number of variants are read and de-duplicated", {
    expect_equal(length(test_profile_1), 375)
    expect_equal(length(test_profile_2), 374)
    expect_equal(length(mcols(test_profile_1)), 17)
    expect_equal(length(mcols(test_profile_2)), 17)
})

test_that("non-standard chromosomes are removed correcty", {
    expect_identical(levels(seqnames(test_profile_1)), c("1", "12"))
    expect_identical(levels(seqnames(test_profile_2)), c("1", "12"))
})
