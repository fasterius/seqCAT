library("seqCAT")
context("Read COSMIC data")

# Input file
cell_lines_file <- system.file("extdata", package = "seqCAT",
                               "subset_CosmicCLP_MutantExport.tsv.gz")
cosmic_file <- system.file("extdata", package = "seqCAT",
    "subset_CosmicCompleteTargetedScreensMutantExport.tsv.gz")

# Read COSMIC data
hct116 <- suppressMessages(read_cosmic(cell_lines_file, sample_name = "HCT116"))
hct116_list <- suppressMessages(list_cosmic(cell_lines_file))
cosmic <- suppressMessages(read_cosmic(cosmic_file, primary_site = "liver"))
cosmic_list <- suppressMessages(list_cosmic(cosmic_file))

# Tests
test_that("the returned object is a dataframe", {
    expect_identical(class(hct116)[1], "data.frame")
    expect_identical(class(cosmic)[1], "data.frame")
})

test_that("the returned number of variants and metadata columns are correct", {
    expect_equal(nrow(hct116), 1)
    expect_equal(nrow(cosmic), 3)
    expect_equal(ncol(hct116[hct116$sample == "HCT116", ]), 41)
    expect_equal(ncol(cosmic[cosmic$sample == "cosmic", ]), 38)
})

test_that("a character list of cell lines are correctly returned", {
    expect_identical(class(hct116_list), "character")
    expect_identical(class(cosmic_list), "character")
})

test_that("correct number of cell lines are listed", {
    expect_equal(length(hct116_list), 164)
    expect_equal(length(cosmic_list), 189)
})

test_that("missing sample names and primary sites are handled correctly", {
    expect_error(read_cosmic(cell_lines_file, sample_name = "HCT111"),
                 "not present in the data")
    expect_error(read_cosmic(cosmic_file, primary_site = "lever"),
                 "not present in the data")
})
