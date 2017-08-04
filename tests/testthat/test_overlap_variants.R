library("CellAuthentication")
context("Overlap variant sets")

# Get GRanges objects containing the variants
data(granges_1)
data(granges_2)

# Create dummy data for testing zero overlaps and missing genotypes
granges_3 <- granges_1[2, ]
granges_3$A1 <- NA

# Tests
test_that("correct number of variants overlap", {
    expect_equal(nrow(overlap_variants(granges_1, granges_2)), 51)
})

zero_1 <- overlap_variants(granges_2, granges_3)
test_that("zero overlaps are handled correctly", {
    expect_equal(length(zero_1[is.na(zero_1)]), 37)
    expect_equal(length(grep("sample", t(zero_1), value = TRUE)), 2)
})

zero_2 <- overlap_variants(granges_1, granges_3)
test_that("missing genotypes are handled correctly", {
    expect_equal(length(zero_2[is.na(zero_2)]), 37)
    expect_equal(length(grep("sample", t(zero_2), value = TRUE)), 2)
})
