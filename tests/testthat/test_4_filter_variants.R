library("CellAuthentication")
context("Optional filtration of overlapping variants")

# Get overlapping variants
data(test_overlaps)

# Tests
test_that("correct number of variants are filtered", {
    expect_equal(nrow(filter_variants(test_overlaps, 15)), 43)
    expect_equal(nrow(filter_variants(test_overlaps, 50)), 24)
})
