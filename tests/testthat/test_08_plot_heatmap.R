library("seqCAT")
context("Plot heatmap for multiple comparisons")

# Load similarities
data(test_similarities)

# Create plot object
plot <- plot_heatmap(test_similarities)
plot_obj <- ggplot2::ggplot_build(plot)

# Get plot parameters based on ggplot2 version
if (utils::packageVersion("ggplot2") < "2.2.1.9000") {
    plot_params <- plot_obj$layout$panel_ranges[[1]]
} else {
    plot_params <- plot_obj$layout$panel_params[[1]]
}

# Tests
test_that("mirrored data is added correctly", {
    expect_equal(plot_obj$data[[2]]$label[1], 91.7)
    expect_equal(plot_obj$data[[2]]$label[9], 37.5)
})

test_that("samples are clustered correctly", {
    expect_identical(plot_params$x.labels[1], "sample3")
    expect_identical(plot_params$x.labels[2], "sample1")
    expect_identical(plot_params$x.labels[3], "sample2")
})
