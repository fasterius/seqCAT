library("seqCAT")
context("Plot impact distributions")

# Load data
data(test_comparison)

# Create plot
plot <- plot_impacts(test_comparison)
plot_obj <- ggplot2::ggplot_build(plot)

# Tests
test_that("correct number of variants per impact are plotted", {
    expect_equal(plot_obj$data[[2]]$label[1], 0)
    expect_equal(plot_obj$data[[2]]$label[4], 1)
    expect_equal(plot_obj$data[[2]]$label[7], 1)
    expect_equal(plot_obj$data[[2]]$label[8], 49)
})

test_that("correct impact distributions are calculated", {
    expect_identical(plot_obj$data[[3]]$label[1], "100.0 %")
    expect_identical(plot_obj$data[[3]]$label[2], " 98.0 %")
    expect_identical(plot_obj$data[[4]]$label[1], "  0.0 %")
    expect_identical(plot_obj$data[[4]]$label[4], "  2.0 %")
    expect_identical(plot_obj$layout$panel_ranges[[1]]$x.labels[4], "MODIFIER")
})
