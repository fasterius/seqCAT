library("CellAuthentication")
context("Overlap SNV profiles")

# Get GRanges objects containing the variants
data(profile_1)
data(profile_2)

# Create dummy data for testing zero overlaps and missing genotypes
profile_3 <- profile_1[2, ]
profile_3$A1 <- NA

# Tests
test_that("correct number of variants overlap", {
    expect_equal(nrow(overlap_profiles(profile_1, profile_2)), 51)
})

zero_1 <- overlap_profiles(profile_2, profile_3)
test_that("profiles with zero overlaps are handled correctly", {
    expect_equal(length(zero_1[is.na(zero_1)]), 37)
    expect_equal(length(grep("sample", t(zero_1), value = TRUE)), 2)
})

zero_2 <- overlap_profiles(profile_1, profile_3)
test_that("missing genotypes are handled correctly", {
    expect_equal(length(zero_2[is.na(zero_2)]), 37)
    expect_equal(length(grep("sample", t(zero_2), value = TRUE)), 2)
})
