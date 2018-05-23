library("seqCAT")
context("Plot lists of known variants")

# Load test data
data(test_variant_list)

# Create plot object
plot <- plot_variant_list(test_variant_list)
plot_obj <- ggplot2::ggplot_build(plot)

# Get plot parameters based on ggplot2 version
if (utils::packageVersion("ggplot2") < "2.2.1.9000") {
    plot_params <- plot_obj$layout$panel_ranges[[1]]
} else {
    plot_params <- plot_obj$layout$panel_params[[1]]
}

# Tests
test_that("correct variant/gene combinations are being plotted", {
    expect_equal(plot_params$x.labels[1], "1:16229 (DDX11L1)")
    expect_equal(plot_params$x.labels[2], "1:16298 (DDX11L1)")
    expect_equal(plot_params$y.labels[1], "sample1")
    expect_equal(plot_params$y.labels[2], "sample2")
})

test_that("genotypes are separated and coloured correctly", {
    expect_equal(plot_obj$data[[1]][1, "fill"], "#999999")
    expect_equal(plot_obj$data[[2]][3, "fill"], "#4e8ce4")
    expect_equal(plot_obj$data[[3]][12, "fill"], "#a6c6f2")
})
