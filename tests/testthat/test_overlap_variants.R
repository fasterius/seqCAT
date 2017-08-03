library("CellAuthentication")
context("Overlap variant sets")

# Files
file1 = system.file("extdata", "variants.sample1.txt",
                   package="CellAuthentication")
file2 = system.file("extdata", "variants.sample2.txt",
                   package="CellAuthentication")

# Read data from files
data1 = read.table(file1, sep="\t", header=TRUE, stringsAsFactors=FALSE)
data2 = read.table(file2, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Create dummy data for testing zero overlaps and missing genotypes
data3 = data1[2, ]
data3$A1 = NA

# Function for converting data frame to GRanges object
convert.df = function(df) {
    gr = GenomicRanges::makeGRangesFromDataFrame(df,
        keep.extra.columns=TRUE, 
        ignore.strand=TRUE,
        seqinfo=NULL, 
        seqnames.field="seqnames",
        start.field="start",
        end.field="end",
        starts.in.df.are.0based=FALSE)
    return(gr)
}

# Convert to GRanges
variants1 = convert.df(data1)
variants2 = convert.df(data2)
variants3 = convert.df(data3)

# Tests
test_that("correct number of variants overlap", {
    expect_equal(nrow(overlap_variants(variants1, variants2)), 51)
})

zero_1 = overlap_variants(variants2, variants3)
test_that("zero overlaps are handled correctly", {
    expect_equal(length(zero_1[is.na(zero_1)]), 37)
    expect_equal(length(grep("sample", t(zero_1), value=TRUE)), 2)
})

zero_2 = overlap_variants(variants1, variants3)
test_that("missing genotypes are handled correctly", {
    expect_equal(length(zero_2[is.na(zero_2)]), 37)
    expect_equal(length(grep("sample", t(zero_2), value=TRUE)), 2)
})
