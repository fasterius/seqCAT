library("CellAuthentication")
context("Read extracted variants")

# Read variants
file_1 <- system.file("extdata",
                      "extract.sample1.txt",
                      package = "CellAuthentication")
file_2 <- system.file("extdata",
                      "extract.sample2.txt",
                      package = "CellAuthentication")

data_1 <- as.data.frame(read_variants(file        = file_1,
                                      sample_name = "sample1"))

data_2 <- as.data.frame(read_variants(file        = file_2,
                                      sample_name = "sample2"))

# Tests
test_that("correct number of variants are read and de-duplicated", {
    expect_equal(nrow(unique(data_1[c("seqnames", "start", "ENSGID")])), 383)
    expect_equal(nrow(unique(data_2[c("seqnames", "start", "ENSGID")])), 382)
})

test_that("only chromosome 1-22, X and Y are read", {
    expect_equal(nrow(data_1[!(data_1$seqnames %in% c(1:22, "X", "Y")), ]), 0)
    expect_equal(nrow(data_2[!(data_2$seqnames %in% c(1:22, "X", "Y")), ]), 0)
})
