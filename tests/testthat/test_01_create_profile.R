library("seqCAT")
context("Creation of SNV profiles")

# Paths and directories
file1 <- system.file("extdata", "test.vcf.gz", package = "seqCAT")
file2 <- system.file("extdata", "test.unannotated.vcf.gz", package = "seqCAT")
file3 <- system.file("extdata", "test.gvcf.gz", package = "seqCAT")
vcf_dir <- system.file("extdata", package = "seqCAT")

# Create individual profiles
suppressMessages(create_profile(vcf_file     = file1,
                                sample       = "sample1",
                                output_file  = "profile_1.txt",
                                filter_depth = 10,
                                filter_vc    = TRUE,
                                python       = FALSE))

suppressMessages(create_profile(vcf_file     = file1,
                                sample       = "sample2",
                                output_file  = "profile_2.txt",
                                filter_depth = 10,
                                filter_vc    = TRUE,
                                python       = FALSE))

suppressMessages(create_profile(vcf_file     = file2,
                                sample       = "sample3",
                                output_file  = "profile_3.txt",
                                filter_depth = 10,
                                filter_vc    = TRUE,
                                python       = FALSE))

# Create profiles in directory
suppressMessages(create_profiles(vcf_dir      = vcf_dir,
                                 output_dir   = ".",
                                 pattern      = "sample1",
                                 recursive    = FALSE,
                                 filter_depth = 10,
                                 filter_vc    = TRUE,
                                 python       = FALSE))

# Read profiles
profile_1 <- read.table(file             = "profile_1.txt",
                        sep              = "\t",
                        header           = TRUE,
                        stringsAsFactors = FALSE)

profile_2 <- read.table(file             = "profile_2.txt",
                        sep              = "\t",
                        header           = TRUE,
                        stringsAsFactors = FALSE)

profile_3 <- read.table(file             = "profile_3.txt",
                        sep              = "\t",
                        header           = TRUE,
                        stringsAsFactors = FALSE)

profile_dir <- read.table(file             = "sample1.profile.txt",
                          sep              = "\t",
                          header           = TRUE,
                          stringsAsFactors = FALSE)

# Tests
test_that("create_profile yields correct dimensions", {
    expect_equal(dim(profile_1), c(428, 18))
    expect_equal(dim(profile_2), c(428, 18))
    expect_equal(dim(profile_dir), c(428, 18))
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
    expect_equal(nrow(profile_1[profile_1$impact == "MODERATE", ]), 4)
    expect_equal(nrow(profile_2[profile_2$impact == "MODERATE", ]), 4)
    expect_equal(nrow(profile_dir[profile_dir$impact == "MODERATE", ]), 4)
    expect_equal(nrow(profile_1[profile_1$impact == "LOW", ]), 0)
    expect_equal(nrow(profile_2[profile_2$impact == "LOW", ]), 0)
    expect_equal(nrow(profile_dir[profile_dir$impact == "LOW", ]), 0)
    expect_equal(nrow(profile_1[profile_1$impact == "MODIFIER", ]), 423)
    expect_equal(nrow(profile_2[profile_2$impact == "MODIFIER", ]), 424)
    expect_equal(nrow(profile_dir[profile_dir$impact == "MODIFIER", ]), 423)
})

test_that("empty calls from multi-sample VCFs are handled correctly", {
    expect_equal(nrow(profile_1[profile_1$rsID == "rs182017058", ]), 0)
    expect_equal(nrow(profile_2[profile_2$rsID == "rs182017058", ]), 3)
})

test_that("VCFs without variant annotations are handled correctly", {
    expect_equal(nrow(profile_3), 99)
})

test_that("VCFs with <NON_REF> alleles are handled properly", {
    expect_warning(create_profile(vcf_file     = file3,
                                  sample       = "sample4",
                                  output_file  = "profile_4.txt",
                                  filter_depth = 10,
                                  filter_vc    = FALSE,
                                  python       = FALSE),
                   "<NON_REF> alleles; input may be a gVCF")
})

# Remove files
file.remove("profile_1.txt")
file.remove("profile_2.txt")
file.remove("profile_3.txt")
file.remove("profile_4.txt")
file.remove("sample1.profile.txt")
