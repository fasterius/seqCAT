library("seqCAT")
context("Listing known variants in SNV profiles")

# Load test profile
data(test_profile_1)
data(test_profile_2)

# Create test variants
test_variants <- data.frame(chr              = c(1, 1),
                            pos              = c(16229, 99999),
                            gene             = c("DDX11L1", "mock_gene"),
                            stringsAsFactors = FALSE)

# List test variants
profiles <- list(test_profile_1, test_profile_2)
variants <- list_variants(profiles, test_variants)

# Tests
test_that("errors for malformed input data are raised correctly", {
    expect_error(list_variants(test_profile_1, "a string"), "not a dataframe")
    expect_error(list_variants(test_profile_1, data.frame(chr = 1)), "'pos'")
    expect_error(list_variants(test_profile_1, data.frame(pos = 1)), "'chr'")
})

test_that("known variants are listed correctly", {
    expect_equal(nrow(variants), 2)
    expect_equal(variants[variants$sample1 == "CA", "gene"], "DDX11L1")
    expect_equal(variants[variants$sample2 == "AA", "gene"], "DDX11L1")
    expect_equal(variants[variants$sample1 == "0", "gene"], "mock_gene")
})
