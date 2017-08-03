library("CellAuthentication")
context("Optional filtration of overlapping variants")

# Get variants
file <- system.file("extdata",
                    "overlaps.txt",
                    package = "CellAuthentication")
variants <- read.table(file,
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = FALSE)

# Tests
test_that("correct number of variants are filtered", {
    expect_equal(nrow(filter_variants(variants, 15)), 43)
    expect_equal(nrow(filter_variants(variants, 50)), 24)
})
