library("CellAuthentication")
context("Plot impact distributions")

# Load data
data(test_comparison)

# Create plot
plot <- plot_impacts(test_comparison)
plot_obj <- ggplot2::ggplot_build(plot)

# Tests
test_that("correct impact distributions are calculated", {
    expect_equal(plot_obj$data[[2]]$label[1], 1)
    expect_equal(plot_obj$data[[2]]$label[2], 50)
    expect_identical(plot_obj$data[[3]]$label[1], "100.0 %")
    expect_identical(plot_obj$layout$panel_ranges[[1]]$x.labels, "MODIFIER")
})
