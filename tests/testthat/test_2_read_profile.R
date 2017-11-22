library("seqCAT")
context("Read SNV profiles")

# Paths
file_1 <- system.file("extdata", "test_profile_1.txt.gz", package = "seqCAT")
file_2 <- system.file("extdata", "test_profile_2.txt.gz", package = "seqCAT")
file_3 <- system.file("extdata", "test_profile_3.txt.gz", package = "seqCAT")

# Read profiles
test_profile_1 <- suppressMessages(read_profile(file        = file_1,
                                                sample_name = "sample1"))

test_profile_2 <- suppressMessages(read_profile(file        = file_2,
                                                sample_name = "sample2"))

test_profile_3 <- suppressMessages(read_profile(file        = file_3,
                                                sample_name = "sample2"))

# Tests
test_that("a GRanges object is returned", {
    expect_identical(class(test_profile_1)[1], "GRanges")
})

test_that("correct number of variants are read and de-duplicated", {
    expect_equal(length(test_profile_1), 376)
    expect_equal(length(test_profile_2), 374)
    expect_equal(length(mcols(test_profile_1)), 17)
    expect_equal(length(mcols(test_profile_2)), 17)
})

test_that("non-standard chromosomes are removed correcty", {
    expect_identical(levels(seqnames(test_profile_1)), c("1", "12"))
    expect_identical(levels(seqnames(test_profile_2)), c("1", "12"))
})

test_that("profiles without variant annotations are handled correctly", {
    expect_equal(length(test_profile_3), 99)
})
