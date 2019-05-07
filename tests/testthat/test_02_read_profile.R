library("seqCAT")
context("Read SNV profiles")

# Paths
file_1 <- system.file("extdata", "test_1.profile.txt.gz", package = "seqCAT")
file_2 <- system.file("extdata", "test_2.profile.txt.gz", package = "seqCAT")
file_3 <- system.file("extdata", "test_3.profile.txt.gz", package = "seqCAT")
file_4 <- system.file("extdata", "test_4.profile.txt.gz", package = "seqCAT")
profile_dir <- system.file("extdata", package = "seqCAT")

# Read profiles
test_1 <- suppressMessages(read_profile(file = file_1))
test_2 <- suppressMessages(read_profile(file = file_2))
test_3 <- suppressMessages(read_profile(file = file_3))
test_4 <- suppressMessages(read_profile(file = file_4))
test_5 <- suppressMessages(read_profile(file = file_4))
test_dir <- suppressMessages(read_profiles(profile_dir = profile_dir))

# Tests
test_that("correct object types are returned", {
    expect_identical(class(test_1)[1], "data.frame")
    expect_identical(class(test_dir)[1], "list")
})

test_that("correct number of variants are read", {
    expect_equal(dim(test_1), c(377, 19))
    expect_equal(dim(test_2), c(375, 19))
    expect_equal(dim(test_dir[[1]]), c(377, 19))
})

test_that("profiles without variant annotations are handled correctly", {
    expect_equal(nrow(test_3), 99)
})
#
# test_that("empty profiles are handled correcty", {
    # expect_equal(nrow(test_4), 1)
    # expect_identical(test_4$sample, "sample4")
# })
