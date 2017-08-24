library("CellAuthentication")
context("Compare many profiles")

# Load profiles
data(test_profile_1)
data(test_profile_2)

# List profiles
profiles <- list(test_profile_1, test_profile_2)

# Comparisons
many <- compare_many(profiles)
one <- compare_many(profiles, profiles[[1]])

# Tests
test_that("the returned object is a list", {
    expect_identical(class(many), "list")
})

test_that("correct number of comparisons are performed", {
    expect_equal(nrow(many[[1]]), 3)
    expect_equal(nrow(one[[1]]), 2)
})

test_that("correct correct summary statistics are calculated", {
    expect_equal(many[[1]][3, 6], 91.2)
    expect_equal(one[[1]][2, 5], 98.0)
})

test_that("comparisons are correctly ordered in the resulting lists", {
    expect_identical(unique(many[[2]][[3]]$sample_1), "sample2")
    expect_identical(unique(one[[2]][[2]]$sample_1), "sample1")
})
