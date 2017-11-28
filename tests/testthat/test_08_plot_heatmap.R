library("seqCAT")
context("Plot heatmap for multiple comparisons")

# Load similarities
data(test_similarities)

# Plots
plot <- plot_heatmap(test_similarities)
plot_obj <- ggplot2::ggplot_build(plot)

# Tests
test_that("mirrored data is added correctly", {
    expect_equal(plot_obj$data[[2]]$label[1], 99.1)
    expect_equal(plot_obj$data[[2]]$label[9], 67.9)
})

test_that("samples are clustered correctly", {
    expect_identical(plot_obj$layout$panel_ranges[[1]]$x.labels[1], "sample_z")
    expect_identical(plot_obj$layout$panel_ranges[[1]]$x.labels[2], "sample_x")
    expect_identical(plot_obj$layout$panel_ranges[[1]]$x.labels[3], "sample_y")
})
