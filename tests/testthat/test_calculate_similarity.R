library("CellAuthentication")
context("Calculation of summary statistics for overlaps")

# Get compared, overlapping variants
data(overlaps)
overlaps <- compare_overlaps(overlaps)

# Tests
test_that("dataframe checking works correctly", {
    expect_error(calculate_similarity(overlaps, 1), "not a dataframe")
})

dummy <- data.frame(sample_1 = 1, sample_2 = 2, dummy = 3)
test_that("dataframe structure checks are working correctly", {
    expect_error(calculate_similarity(overlaps, dummy), "correct structure")
})

similarity <- calculate_similarity(overlaps)
test_that("correct statistics are calculated", {
    expect_equal(similarity[1, "overlaps"], 51)
    expect_equal(similarity[1, "matches"], 50)
    expect_equal(similarity[1, "concordance"], 98.0)
    expect_equal(similarity[1, "similarity_score"], 89.5)
})

similarity <- calculate_similarity(overlaps, similarity)
test_that("iterative runs are correctly added to each other", {
    expect_equal(similarity[2, "similarity_score"], 89.5)
})
