library("CellAuthentication")
context("VCF variant extraction")

# Files
vcf_file = system.file("extdata", "example.vcf", package="CellAuthentication")
extract = 'extract.txt'

# Extract variants (sample 1)
suppressMessages(extract_variants(vcf_file, "sample1", extract,
                                  filter_depth=10, python=FALSE))
variants_1 = read.table(extract, sep="\t", header=TRUE,
                      stringsAsFactors=FALSE)

# Extract variants (sample 2)
suppressMessages(extract_variants(vcf_file, "sample2", extract,
                                  filter_depth=10, python=FALSE))
variants_2 = read.table(extract, sep="\t", header=TRUE,
                      stringsAsFactors=FALSE)

# Remove extract file
file.remove(extract)

# Tests
test_that("extract_variants yields correct number of columns", {
    expect_equal(ncol(variants_1), 18)
    expect_equal(ncol(variants_2), 18)
})

test_that("extract_variants yields correct number of variants", {
    expect_equal(nrow(variants_1), 433)
    expect_equal(nrow(variants_2), 431)
})

test_that("only variants passing the depth threshold are extracted", {
    expect_equal(nrow(variants_1[variants_1$DP == 0, ]), 0)
    expect_equal(nrow(variants_2[variants_2$DP == 0, ]), 0)
})

test_that("no indels are extracted", {
    expect_equal(nrow(variants_1[nchar(variants_1$REF) != 1, ]), 0)
    expect_equal(nrow(variants_1[nchar(variants_1$ALT) != 1, ]), 0)
    expect_equal(nrow(variants_2[nchar(variants_2$REF) != 1, ]), 0)
    expect_equal(nrow(variants_2[nchar(variants_2$ALT) != 1, ]), 0)
})

test_that("the MT chromosome only exists in sample 1", {
    expect_equal(nrow(variants_1[variants_1$chr == "MT", ]), 1)
    expect_equal(nrow(variants_2[variants_2$chr == "MT", ]), 0)
})
