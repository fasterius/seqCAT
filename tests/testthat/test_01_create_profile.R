library("seqCAT")
context("Creation of SNV profiles")

# Paths and directories
file1 <- system.file("extdata", "test.vcf.gz", package = "seqCAT")
file2 <- system.file("extdata", "test.unannotated.vcf.gz", package = "seqCAT")
file3 <- system.file("extdata", "test.gvcf.gz", package = "seqCAT")
vcf_dir <- system.file("extdata", package = "seqCAT")

# Create individual profiles
profile_1 <- suppressMessages(create_profile(vcf_file  = file1,
                                             sample    = "sample1",
                                             min_depth = 10,
                                             filter_vc = TRUE,
                                             filter_mt = FALSE,
                                             filter_ns = TRUE,
                                             filter_gd = TRUE,
                                             filter_pd = FALSE))

profile_2 <- suppressMessages(create_profile(vcf_file  = file1,
                                             sample    = "sample2",
                                             min_depth = 10,
                                             filter_vc = TRUE,
                                             filter_mt = FALSE,
                                             filter_ns = FALSE,
                                             filter_gd = FALSE,
                                             filter_pd = TRUE))

profile_3 <- suppressMessages(create_profile(vcf_file  = file2,
                                             sample    = "sample3",
                                             min_depth = 10,
                                             filter_vc = TRUE,
                                             filter_mt = FALSE,
                                             filter_ns = TRUE,
                                             filter_gd = FALSE,
                                             filter_pd = TRUE))

# Create profiles in directory
profile_dir <- suppressMessages(create_profiles(vcf_dir   = vcf_dir,
                                                min_depth = 10,
                                                filter_vc = TRUE,
                                                filter_mt = FALSE,
                                                filter_ns = TRUE,
                                                filter_gd = TRUE,
                                                filter_pd = FALSE,
                                                pattern   = "sample1",
                                                recursive = FALSE))[[1]]

# Tests
test_that("create_profile yields correct dimensions", {
    expect_equal(dim(profile_1), c(374, 20))
    expect_equal(dim(profile_2), c(52, 20))
    expect_equal(dim(profile_dir), c(374, 20))
})

test_that("the correct variants across impact categories are extracted", {
    expect_equal(nrow(profile_1[profile_1$impact == "HIGH", ]), 1)
    expect_equal(nrow(profile_2[profile_2$impact == "HIGH", ]), 0)
    expect_equal(nrow(profile_dir[profile_dir$impact == "HIGH", ]), 1)
    expect_equal(nrow(profile_1[profile_1$impact == "MODERATE", ]), 1)
    expect_equal(nrow(profile_2[profile_2$impact == "MODERATE", ]), 1)
    expect_equal(nrow(profile_dir[profile_dir$impact == "MODERATE", ]), 1)
    expect_equal(nrow(profile_1[profile_1$impact == "LOW", ]), 0)
    expect_equal(nrow(profile_2[profile_2$impact == "LOW", ]), 0)
    expect_equal(nrow(profile_dir[profile_dir$impact == "LOW", ]), 0)
    expect_equal(nrow(profile_1[profile_1$impact == "MODIFIER", ]), 372)
    expect_equal(nrow(profile_2[profile_2$impact == "MODIFIER", ]), 51)
    expect_equal(nrow(profile_dir[profile_dir$impact == "MODIFIER", ]), 372)
})

test_that("empty calls from multi-sample VCFs are handled correctly", {
    expect_equal(nrow(profile_1[profile_1$rsID == "rs182017058", ]), 0)
    expect_equal(nrow(profile_2[profile_2$rsID == "rs182017058", ]), 1)
})

test_that("VCFs without variant annotations are handled correctly", {
    expect_equal(nrow(profile_3), 99)
})

test_that("Samples not present in the VCF file are handled properly", {
    expect_error(create_profile(vcf_file     = file1,
                                sample       = "sampleX",
                                min_depth    = 10,
                                filter_vc    = TRUE),
                 "is not present in the VCF file")
})
