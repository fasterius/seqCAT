library("seqCAT")
context("Write SNV profiles")

# Test data
data(test_profile_1)

# Tests
test_that("unsupported format specifications are handled correctly", {
    expect_error(write_profile(test_profile_1,
                               "test.bef"),
                 "Unsupported format specification \"")
})
