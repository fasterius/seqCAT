library("CellAuthentication")
context("Overlap variant sets")

# Get variants
data(variants_1)
data(variants_2)
# file1 <- system.file("extdata",
                     # "variants.sample1.txt",
                     # package = "CellAuthentication")
# file2 <- system.file("extdata",
                     # "variants.sample2.txt",
                     # package = "CellAuthentication")
#
# # Read data from files
# data1 <- read.table(file             = file1,
                    # sep              = "\t",
                    # header           = TRUE,
                    # stringsAsFactors = FALSE)
#
# data2 <- read.table(file             = file2,
                    # sep              = "\t",
                    # header           = TRUE,
                    # stringsAsFactors = FALSE)

# Create dummy data for testing zero overlaps and missing genotypes
variants_3 <- variants_1[2, ]
variants_3$A1 <- NA

# Function for converting data frame to GRanges object
convert_df <- function(variants) {
    gr <- GenomicRanges::makeGRangesFromDataFrame(variants,
        keep.extra.columns      = TRUE,
        ignore.strand           = TRUE,
        seqinfo                 = NULL,
        seqnames.field          = "seqnames",
        start.field             = "start",
        end.field               = "end",
        starts.in.df.are.0based = FALSE)
    return(gr)
}

# Convert to GRanges
variants_1 <- convert_df(variants_1)
variants_2 <- convert_df(variants_2)
variants_3 <- convert_df(variants_3)

# Tests
test_that("correct number of variants overlap", {
    expect_equal(nrow(overlap_variants(variants_1, variants_2)), 51)
})

zero_1 <- overlap_variants(variants_2, variants_3)
test_that("zero overlaps are handled correctly", {
    expect_equal(length(zero_1[is.na(zero_1)]), 37)
    expect_equal(length(grep("sample", t(zero_1), value = TRUE)), 2)
})

zero_2 <- overlap_variants(variants_1, variants_3)
test_that("missing genotypes are handled correctly", {
    expect_equal(length(zero_2[is.na(zero_2)]), 37)
    expect_equal(length(grep("sample", t(zero_2), value = TRUE)), 2)
})
