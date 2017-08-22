library("CellAuthentication")
context("Overlap and compare SNV profiles")

# Get GRanges objects containing the variants
data(test_profile_1)
data(test_profile_2)

# Create dummy data for testing zero overlaps and missing genotypes
test_profile_3 <- test_profile_1[2, ]
test_profile_3$A1 <- NA

# Comparisons
compared <- compare_profiles(test_profile_1, test_profile_2)
zero_1 <- compare_profiles(test_profile_2, test_profile_3)
zero_2 <- compare_profiles(test_profile_1, test_profile_3)

# Tests
test_that("correct number of variants overlap", {
    expect_equal(nrow(compared), 51)
})

test_that("correct number of matches and mismatches are found", {
    expect_equal(nrow(compared), 51)
    expect_equal(nrow(compared[compared$match == "match", ]), 50)
    expect_equal(nrow(compared[compared$match == "mismatch", ]), 1)
})

alleles <- paste(c("A1", "A1", "A2", "A2"),
                 c("sample_1", "sample_2"),
                 sep = ".")
test_that("no missing genotypes are found", {
    expect_equal(nrow(compared[complete.cases(compared[, alleles]), ]), 51)
})

test_that("profiles with zero overlaps are handled correctly", {
    expect_equal(length(zero_1[is.na(zero_1)]), 15)
    expect_equal(length(grep("sample", t(zero_1), value = TRUE)), 2)
})

test_that("missing genotypes are handled correctly", {
    expect_equal(length(zero_2[is.na(zero_2)]), 15)
    expect_equal(length(grep("sample", t(zero_2), value = TRUE)), 2)
})
