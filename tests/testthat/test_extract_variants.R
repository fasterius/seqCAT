library("CellAuthentication")
context("VCF variant extraction")

# Files
file <- system.file("extdata",
                    "example.vcf",
                    package = "CellAuthentication")

# Extract variants
suppressMessages(extract_variants(vcf_file     = file,
                                  sample       = "sample1",
                                  output_file  = "extract_1.txt",
                                  filter_depth = 10,
                                  python       = FALSE))

suppressMessages(extract_variants(vcf_file     = file,
                                  sample       = "sample2",
                                  output_file  = "extract_2.txt",
                                  filter_depth = 10,
                                  python       = FALSE))

# Read files
variants_1 <- read.table(file             = "extract_1.txt",
                         sep              = "\t",
                         header           = TRUE,
                         stringsAsFactors = FALSE)

variants_2 <- read.table(file             = "extract_2.txt",
                         sep              = "\t",
                         header           = TRUE,
                         stringsAsFactors = FALSE)

# Remove extract file
file.remove("extract_1.txt")
file.remove("extract_2.txt")

# Tests
test_that("extract_variants yields correct dimensions", {
    expect_equal(dim(variants_1), c(432, 18))
    expect_equal(dim(variants_2), c(430, 18))
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
