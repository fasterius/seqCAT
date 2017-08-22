library("CellAuthentication")
context("Optional filtration of overlapping variants")

# Get overlapping variants
data(test_comparison)

# Tests
test_that("correct number of variants are filtered", {
    expect_equal(nrow(filter_variants(test_comparison, 15)), 43)
    expect_equal(nrow(filter_variants(test_comparison, 50)), 24)
})
