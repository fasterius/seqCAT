library("CellAuthentication")
context("Compare SNV profile overlaps")

# Get overlapping variants
data(overlaps)

# Tests
compared <- compare_overlaps(overlaps)
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
