library("CellAuthentication")
context("Read extracted variants")

# Files
file1 = system.file("extdata", "extract.sample1.txt", 
                    package="CellAuthentication")
file2 = system.file("extdata", "extract.sample2.txt", 
                    package="CellAuthentication")
file3 = system.file("extdata", "variants.sample1.txt",
                    package="CellAuthentication")
file4 = system.file("extdata", "variants.sample2.txt",
                    package="CellAuthentication")

# Read variants
data1 = as.data.frame(read_variants(file1, "sample1"))
data2 = as.data.frame(read_variants(file2, "sample2"))
expected1 = read.table(file3, sep="\t", header=TRUE, stringsAsFactors=FALSE)
expected2 = read.table(file4, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Tests
test_that("variants are correctly de-duplicated", {
    expect_equal(nrow(unique(data1[c("seqnames", "start", "ENSGID")])), 383)
    expect_equal(nrow(unique(data2[c("seqnames", "start", "ENSGID")])), 382)
})

test_that("only chromosome 1-22, X and Y are read", {
    expect_equal(nrow(data1[!(data1$seqnames %in% c(1:22, "X", "Y")), ]), 0)
    expect_equal(nrow(data2[!(data2$seqnames %in% c(1:22, "X", "Y")), ]), 0)
})
