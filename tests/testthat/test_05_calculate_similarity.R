library("seqCAT")
context("Calculation of summary statistics for overlaps")

# Get compared, overlapping variants
data(test_comparison)

# Calculate similarities
similarity <- calculate_similarity(test_comparison)
similarity <- calculate_similarity(test_comparison, similarity)

# Dummy for structure test
dummy <- data.frame(sample_1 = 1, sample_2 = 2, dummy = 3)

# Tests
test_that("dataframe structure checks are working correctly", {
    expect_error(calculate_similarity(test_comparison, similarity = 1),
                 "not a dataframe")
    expect_error(calculate_similarity(test_comparison, dummy), "structure")
})

test_that("correct statistics are calculated", {
    expect_equal(similarity[1, "overlaps"], 51)
    expect_equal(similarity[1, "matches"], 50)
    expect_equal(similarity[1, "concordance"], 98.0)
    expect_equal(similarity[1, "similarity_score"], 89.5)
})

test_that("iterative runs are correctly added to each other", {
    expect_equal(similarity[2, "similarity_score"], 89.5)
})
