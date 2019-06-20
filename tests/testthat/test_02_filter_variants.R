library("seqCAT")
context("Variant filtration")

# Get test profiles
data(test_profile_1)
data(test_profile_2)

# Create additional test profiles
test_profile_3 <- test_profile_1
test_profile_3$ALT <- "<NON_REF>"
test_profile_4 <- test_profile_1
test_profile_4[1, "ALT"] <- "<NON_REF>"
test_profile_5 <- test_profile_1
test_profile_5$FILTER <- "."
test_profile_6 <- test_profile_1
test_profile_6$ENSGID <- NULL

# Tests
test_that("correct number of variants are filtered on depth", {
    expect_equal(nrow(filter_variants(test_profile_1, min_depth = 15)), 328)
})

test_that("the MT chromosome only exists in sample 1", {
    expect_equal(nrow(test_profile_1[test_profile_1$chr == "MT", ]), 1)
    expect_equal(nrow(test_profile_2[test_profile_2$chr == "MT", ]), 0)
})

test_that("no indels are present in the final profiles", {
    expect_equal(nrow(test_profile_1[nchar(test_profile_1$REF) != 1, ]), 0)
    expect_equal(nrow(test_profile_1[nchar(test_profile_1$ALT) != 1, ]), 0)
    expect_equal(nrow(test_profile_2[nchar(test_profile_2$REF) != 1, ]), 0)
    expect_equal(nrow(test_profile_2[nchar(test_profile_2$ALT) != 1, ]), 0)
})

test_that("VCFs with <NON_REF> alleles are handled properly", {
    expect_error(filter_variants(test_profile_3),
                   "<NON_REF> alleles; input may be a gVCF")
    expect_warning(filter_variants(test_profile_4),
                   "<NON_REF> alleles; input may be a gVCF")
})

test_that("VCFs without FILTER data are handled correctly", {
    expect_error(filter_variants(test_profile_5, filter_vc = TRUE),
                 "VCF contains no FILTER data; please")
})

test_that("Profiles with zero variants after filtering are handled properly", {
    expect_error(filter_variants(test_profile_1, min_depth = 10000),
                 "No variants left after filtering with the current")
})

test_that("De-duplication yields the correct number of variants", {
    expect_equal(nrow(filter_duplicates(test_profile_1, filter_pd = TRUE)), 54)
})

test_that("Profiles without gene-level information are handled correctly", {
    expect_error(filter_duplicates(test_profile_6, filter_gd = TRUE),
                 "No ENSGID gene data available")
})
