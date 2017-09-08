library("seqCAT")
context("Read COSMIC data")

# Input file
file <- system.file("extdata",
                    "subset_CosmicCLP_MutantExport.tsv.gz",
                    package = "seqCAT")

# Read COSMIC data
cosmic <- suppressMessages(read_cosmic(file, "HCT116"))
cosmic_list <- suppressMessages(list_cosmic(file))

# Tests
test_that("the returned object is a GRanges object", {
    expect_identical(class(cosmic)[1], "GRanges")
})

test_that("the returned number of variants and metadata columns are correct", {
    expect_equal(length(cosmic), 1)
    expect_equal(length(mcols(cosmic[cosmic$sample == "HCT116", ])), 12)
})

test_that("a character list of cell lines are correctly returned", {
    expect_identical(class(cosmic_list), "character")
})

test_that("correct number of cell lines are listed", {
    expect_equal(length(cosmic_list), 164)
})
