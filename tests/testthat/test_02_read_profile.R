library("seqCAT")
context("Read SNV profiles")

# Paths
file_1 <- system.file("extdata", "test_1.profile.txt.gz", package = "seqCAT")
file_2 <- system.file("extdata", "test_2.profile.txt.gz", package = "seqCAT")
file_3 <- system.file("extdata", "test_3.profile.txt.gz", package = "seqCAT")
file_4 <- system.file("extdata", "test_4.profile.txt.gz", package = "seqCAT")
profile_dir <- system.file("extdata", package = "seqCAT")

# Read profiles
test_1 <- suppressMessages(read_profile(file        = file_1,
                                        sample_name = "sample1"))

test_2 <- suppressMessages(read_profile(file        = file_2,
                                        sample_name = "sample2"))

test_3 <- suppressMessages(read_profile(file        = file_3,
                                        sample_name = "sample3"))

test_4 <- suppressMessages(read_profile(file        = file_4,
                                        sample_name = "sample4"))

test_5 <- suppressMessages(read_profile(file        = file_4,
                                        sample_name = "sample5",
                                        remove_mt   = FALSE))

test_dir <- suppressMessages(read_profiles(profile_dir = profile_dir))

# Tests
test_that("correct object types are returned", {
    expect_identical(class(test_1)[1], "GRanges")
    expect_identical(class(test_dir)[1], "list")
})

test_that("correct number of variants are read and de-duplicated", {
    expect_equal(length(test_1), 376)
    expect_equal(length(test_2), 374)
    expect_equal(length(test_dir[[1]]), 376)
    expect_equal(length(mcols(test_1)), 17)
    expect_equal(length(mcols(test_2)), 17)
    expect_equal(length(mcols(test_dir[[1]])), 17)
})

test_that("non-standard chromosomes are removed correcty", {
    expect_identical(levels(seqnames(test_1)), c("1", "12"))
    expect_identical(levels(seqnames(test_2)), c("1", "12"))
    expect_identical(levels(seqnames(test_dir[[1]])), c("1", "12"))
})

test_that("profiles without variant annotations are handled correctly", {
    expect_equal(length(test_3), 99)
})

test_that("empty profiles are handled correcty", {
    expect_equal(length(test_4), 1)
    expect_identical(test_4$sample, "sample4")
})

test_that("mitochondrial variants are read correctly", {
    expect_equal(length(test_5), 63)
})
