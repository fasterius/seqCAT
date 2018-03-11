#' @title Plot similarity heatmap
#'
#' @description Plot a heatmap of similarities from many-to-many SNV profile
#'  comparisons.
#'
#' @details This function creates publication-ready plots of heatmaps for
#' many-to-many  sample comparisons, taking a long-format dataframe containing
#' the summary statistics of each comparison as input.
#'
#' @export
#' @rdname plot_heatmap
#' @param similarities The long-format dataframe containing the data.
#' @param annotate Annotate each cell with the score (boolean).
#' @param annotate_size Text size of the annotations (numeric).
#' @param legend Show a legend for the colour gradient (boolean).
#' @param legend_size Height and width of the legend (vector of two integers).
#' @param cluster Cluster the samples based on similarity (boolean).
#' @param text_size Text size for axes, labels and legend (numeric).
#' @param limits The limits for the colour gradient (vector of four integers).
#' @param colour The main colour to use for the gradient (character).
#' @return A ggplot2 graphical object.
#'
#' @examples
#' # Load test similarities
#' data(test_similarities)
#'
#' # Plot a similarity heatmap
#' heatmap <- plot_heatmap(test_similarities)
plot_heatmap <- function(similarities,
                         cluster       = TRUE,
                         annotate      = TRUE,
                         annotate_size = 9,
                         legend        = TRUE,
                         legend_size   = c(36, 8),
                         limits        = c(0, 50, 90, 100),
                         text_size     = 14,
                         colour        = "#1954A6") {

    # Add mirrored data
    similarities <- mirror(similarities)

    # Cluster data (if applicable)
    if (cluster) {
        similarities <- cluster_data(similarities)
    }

    # Plot
    heatmap <- ggplot2::ggplot(similarities,
                   ggplot2::aes_string(x    = "sample_1",
                                       y    = "sample_2",
                                       fill = "similarity_score")) +
        ggplot2::geom_tile(colour = "white", size = 0.3) +
        ggplot2::coord_equal() +
        ggplot2::theme(text        = ggplot2::element_text(size = text_size),
                       legend.key.height = grid::unit(legend_size[1], "pt"),
                       legend.key.width  = grid::unit(legend_size[2] *
                                                      ggplot2::.pt, "pt"),
                       axis.ticks        = ggplot2::element_blank(),
                       axis.text.x       = ggplot2::element_text(angle = 90,
                                                                 hjust = 1,
                                                                 vjust = 0.5),
                       panel.background  = ggplot2::element_blank()) +
        ggplot2::labs(x    = NULL,
                      y    = NULL,
                      fill = "Score") +
        ggplot2::scale_fill_gradientn(limits  = c(0, 100),
                                      values  = scales::rescale(limits),
                                      colours = c("#FFFFFF", "#FFFFFF",
                                                  "#808080", colour))

    # Add annotation (if applicable)
    if (annotate) {
        heatmap <- heatmap +
            ggplot2::geom_text(colour = "#FFFFFF",
                               size   = annotate_size / ggplot2::.pt,
                               ggplot2::aes_string(label = "similarity_score"))
    }

    # Add legend (if applicable)
    if (!legend) {
        heatmap <- heatmap +
            ggplot2::theme(legend.position = "none")
    }

    # Return plot
    return(heatmap)

}

# Function for adding mirrored data to complete the heatmap
mirror <- function(similarities) {

    # Get the non-sample column names
    columns <- names(similarities)
    columns <- columns[!(columns %in% c("sample_1", "sample_2"))]

    # Add "y vs. x" data not available in the similarities dataframe
    for (n in seq_len(nrow(similarities))) {

        # Get samples
        sample_1 <- similarities[n, "sample_1"]
        sample_2 <- similarities[n, "sample_2"]

        # Check if reverse exists
        reverse <- similarities[similarities$sample_1 == sample_2 &
                                similarities$sample_2 == sample_1, ]
        if (nrow(reverse) == 0) {

            # Get number of rows
            rows <- nrow(similarities)

            # Add reverse data (sample columns)
            similarities[rows + 1, "sample_1"] <- sample_2
            similarities[rows + 1, "sample_2"] <- sample_1

            # Add reverse data (non-sample columns)
            for (col in columns) {
                similarities[rows + 1, col] <- similarities[n, col]
            }
        }
    }

    # Return with the mirrored data
    return(similarities)
}

# Function for clustering similarity data
cluster_data <- function(similarities) {

    # Get samples
    samples <- unique(similarities$sample_1)
    samples_n <- length(samples)

    # Create empty matrix-like dataframe
    similarities_mat <- as.data.frame(matrix(0,
                                             ncol = samples_n,
                                             nrow = samples_n))
    row.names(similarities_mat) <- samples
    names(similarities_mat) <- samples

    # Add data to dataframe
    for (sample_1 in samples) {
        for (sample_2 in samples) {
            similarities_mat[sample_1, sample_2] <-
                similarities[similarities$sample_1 == sample_1 &
                             similarities$sample_2 == sample_2,
                                 "similarity_score"]
        }
    }

    # Clustering
    clust <- stats::hclust(stats::dist(similarities_mat))
    order <- names(similarities_mat[clust$order])

    # Order data according to clustering
    similarities$sample_1 <- factor(similarities$sample_1, levels = order)
    similarities$sample_2 <- factor(similarities$sample_2, levels = order)

    # Return clustered dataframe
    return(similarities)
}
