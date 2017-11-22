library("seqCAT")
context("Overlap and compare SNV profiles")

# Get GRanges objects containing the variants
data(test_profile_1)
data(test_profile_2)
data(test_profile_3)

# Create dummy data for testing zero overlaps and missing genotypes
test_profile_4 <- test_profile_1[2, ]
test_profile_4$A1 <- NA
test_profile_4$sample <- "sample3"

# Comparisons
compared <- suppressMessages(compare_profiles(test_profile_1, test_profile_2))
ann_one <- suppressMessages(compare_profiles(test_profile_1, test_profile_3))
ann_none <- suppressMessages(compare_profiles(test_profile_3, test_profile_3))
zero_1 <- suppressMessages(compare_profiles(test_profile_2, test_profile_4))
zero_2 <- suppressMessages(compare_profiles(test_profile_1, test_profile_4))

# Tests
test_that("correct number of variants overlap", {
    expect_equal(nrow(compared), 51)
})

test_that("correct number of matches and mismatches are found", {
    expect_equal(nrow(compared), 51)
    expect_equal(nrow(compared[compared$match == "match", ]), 50)
    expect_equal(nrow(compared[compared$match == "mismatch", ]), 1)
})

test_that("comparisons without annotations are handled correctly", {
    expect_equal(nrow(ann_one), 2)
    expect_equal(nrow(ann_none), 99)
})

alleles <- paste(c("A1", "A1", "A2", "A2"),
                 c("sample1", "sample2"),
                 sep = ".")
test_that("no missing genotypes are found", {
    expect_equal(nrow(compared[complete.cases(compared[, alleles]), ]), 51)
})

test_that("profiles with zero overlaps are handled correctly", {
    expect_identical(zero_1[1, "match"], "no overlaps")
    expect_equal(length(grep("sample", t(zero_1), value = TRUE)), 2)
})

test_that("missing genotypes are handled correctly", {
    expect_identical(zero_2[1, "match"], "no overlaps")
    expect_equal(length(grep("sample", t(zero_2), value = TRUE)), 2)
})
