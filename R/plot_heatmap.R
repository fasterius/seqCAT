#' Plot a heatmap of multiple sample comparisons
#'
#' This function creates publication-ready plots of heatmaps for many-to-many
#' sample comparisons, taking a long-format dataframe containing the summary
#' statistics of each comparison as input.
#'
#' @export
#' @rdname plot_heatmap
#' @param similarities The long-format dataframe to be plotted
#' @param annotate Annotate each cell with the score
#' @param annotate_size The size of the annotations
#' @param legend Show a legend for the colour gradient
#' @param cluster Cluster the samples based on similarity
#' @param limits The limits for the colour gradient
#' @param colour The main colour to use for the gradient
#' @return A ggplot2 graphical object
#' @examples
#' data(test_similarities)
#' plot_heatmap(test_similarities)
plot_heatmap <- function(similarities,
                         annotate      = TRUE,
                         annotate_size = 5,
                         legend        = TRUE,
                         cluster       = TRUE,
                         limits        = c(0, 50, 90, 100),
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
        ggplot2::theme(axis.ticks       = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       axis.text.x      = ggplot2::element_text(angle = 90,
                                                                hjust = 1,
                                                                vjust = 0.5)) +
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
                               size   = annotate_size,
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

    # Add "y vs. x" data not available in the similarities dataframe
    for (n in seq(1, nrow(similarities))) {

        # Get samples
        sample_1 <- similarities[n, "sample_1"]
        sample_2 <- similarities[n, "sample_2"]

        # Check if reverse exists
        reverse <- similarities[similarities$sample_1 == sample_2 &
                                similarities$sample_2 == sample_1, ]
        if (nrow(reverse) == 0) {

            # Get number of rows
            rows <- nrow(similarities)

            # Add reverse data
            similarities[rows + 1, "sample_1"] <- sample_2
            similarities[rows + 1, "sample_2"] <- sample_1
            similarities[rows + 1, "overlaps"] <- similarities[n, "overlaps"]
            similarities[rows + 1, "matches"] <- similarities[n, "matches"]
            similarities[rows + 1, "concordance"] <-
                similarities[n, "concordance"]
            similarities[rows + 1, "similarity_score"] <-
                similarities[n, "similarity_score"]

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
