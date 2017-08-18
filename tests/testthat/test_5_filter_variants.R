library("CellAuthentication")
context("Optional filtration of overlapping variants")

# Get overlapping variants
data(overlaps)

# Tests
test_that("correct number of variants are filtered", {
    expect_equal(nrow(filter_variants(overlaps, 15)), 43)
    expect_equal(nrow(filter_variants(overlaps, 50)), 24)
})
