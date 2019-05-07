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
                                             filter    = TRUE,
                                             remove_mt = FALSE))

profile_2 <- suppressMessages(create_profile(vcf_file  = file1,
                                             sample    = "sample2",
                                             min_depth = 10,
                                             filter    = TRUE,
                                             remove_mt = FALSE))

profile_3 <- suppressMessages(create_profile(vcf_file  = file2,
                                             sample    = "sample3",
                                             min_depth = 10,
                                             filter    = TRUE,
                                             remove_mt = FALSE))

# Create profiles in directory
profile_dir <- suppressMessages(create_profiles(vcf_dir   = vcf_dir,
                                                pattern   = "sample1",
                                                recursive = FALSE,
                                                min_depth = 10,
                                                filter    = TRUE,
                                                remove_mt = FALSE))[[1]]

# Tests
test_that("create_profile yields correct dimensions", {
    expect_equal(dim(profile_1), c(377, 19))
    expect_equal(dim(profile_2), c(375, 19))
    expect_equal(dim(profile_dir), c(377, 19))
})

test_that("only variants passing the depth threshold are extracted", {
    expect_equal(nrow(profile_1[profile_1$DP == 0, ]), 0)
    expect_equal(nrow(profile_2[profile_2$DP == 0, ]), 0)
    expect_equal(nrow(profile_dir[profile_dir$DP == 0, ]), 0)
})

test_that("no indels are extracted", {
    expect_equal(nrow(profile_1[nchar(profile_1$REF) != 1, ]), 0)
    expect_equal(nrow(profile_1[nchar(profile_1$ALT) != 1, ]), 0)
    expect_equal(nrow(profile_2[nchar(profile_2$REF) != 1, ]), 0)
    expect_equal(nrow(profile_2[nchar(profile_2$ALT) != 1, ]), 0)
    expect_equal(nrow(profile_dir[nchar(profile_dir$REF) != 1, ]), 0)
    expect_equal(nrow(profile_dir[nchar(profile_dir$ALT) != 1, ]), 0)
})

test_that("the MT chromosome only exists in sample 1", {
    expect_equal(nrow(profile_1[profile_1$chr == "MT", ]), 1)
    expect_equal(nrow(profile_2[profile_2$chr == "MT", ]), 0)
    expect_equal(nrow(profile_dir[profile_dir$chr == "MT", ]), 1)
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
    expect_equal(nrow(profile_1[profile_1$impact == "MODIFIER", ]), 375)
    expect_equal(nrow(profile_2[profile_2$impact == "MODIFIER", ]), 374)
    expect_equal(nrow(profile_dir[profile_dir$impact == "MODIFIER", ]), 375)
})

test_that("empty calls from multi-sample VCFs are handled correctly", {
    expect_equal(nrow(profile_1[profile_1$rsID == "rs182017058", ]), 0)
    expect_equal(nrow(profile_2[profile_2$rsID == "rs182017058", ]), 1)
})

test_that("VCFs without variant annotations are handled correctly", {
    expect_equal(nrow(profile_3), 99)
})

test_that("VCFs without FILTER data are handled correctly", {
    expect_error(create_profile(vcf_file     = file3,
                                sample       = "sample4",
                                min_depth    = 10,
                                filter       = TRUE),
                 "VCF contains no FILTER data; please")
})

test_that("VCFs with <NON_REF> alleles are handled properly", {
    expect_warning(create_profile(vcf_file     = file3,
                                  sample       = "sample4",
                                  min_depth    = 10,
                                  filter       = FALSE),
                   "<NON_REF> alleles; input may be a gVCF")
})
