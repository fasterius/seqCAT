#' @title Plot SNV impact distribution
#'
#' @description \link{plot_impacts} plots SNV impact distributions for a
#'  binary comparison
#'
#' @details This function creates publication-ready plots of the impact
#' distribution from a binary dataset comparison across the matched/mismatched
#' SNVs.
#'
#' @export
#' @rdname plot_impacts
#' @param comparison The profile comparison to be plotted
#' @param annotate Annotate each category
#' @param legend Show/hide legend
#' @param palette Colour palette for filling of bars
#' @return A ggplot2 graphical object
#'
#' @examples
#' # Load test comparison data
#' data(test_comparison)
#'
#' # Plot the impact distribution
#' impacts <- plot_impacts(test_comparison)
plot_impacts <- function(comparison,
                         annotate = TRUE,
                         legend = TRUE,
                         palette = c("#0D2D59", "#1954A6")) {

    # Factorise match and impact columns
    matches <- c("match", "mismatch")
    impacts <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
    comparison$match <- factor(comparison$match, levels = matches)
    comparison$impact <- factor(comparison$impact, levels = impacts)

    # Calculate impact distribution
    data <- dplyr::group_by_(comparison, "match", "impact")
    data <- dplyr::summarise_(data, count = ~n())
    data <- dplyr::mutate_(data, .dots = stats::setNames(
        lazyeval::interp("count / sum(count) * 100"), "prop"))

    # Add zeroes to empty groups
    for (impact in impacts) {
        for (match in matches) {

            # Check current combination for existing data
            current <- data[data$impact == impact &
                            data$match == match, ]$count

            # Add to non-existing data
            if (length(current) == 0) {

                data[nrow(data) + 1, "match"] <- match
                data[nrow(data), "impact"] <- impact
                data[nrow(data), "count"] <- 0
                data[nrow(data), "prop"] <- 0
            }
        }
    }

    # Create plot object
    gg <- ggplot2::ggplot(data, ggplot2::aes_string(x    = "impact",
                                              y    = "prop",
                                              fill = "match")) +
        ggplot2::geom_bar(stat     = "identity",
                          position = "dodge",
                          colour   = "#000000",
                          size     = 0.3) +
        ggplot2::theme_bw() +
        ggplot2::labs(x    = "Impact category",
                      y    = "Proportion of SNVs in category (%)",
                      fill = "") +
        ggplot2::scale_fill_manual(values = palette,
                                   labels = c("Match", "Mismatch")) +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

    # Add text annotation (if applicable)
    if (annotate) {

        # Create labels
        data$percent <- paste0(format(data$prop, nsmall = 1,
                                        digits = 1), " %")

        # Add annotations to plot
        gg <- gg +
            ggplot2::geom_text(
                ggplot2::aes_string(label = "count"),
                position = ggplot2::position_dodge(width = 0.9),
                vjust    = -0.5,
                size     = 2.5,
                colour   = "#000000") +
            ggplot2::geom_text(
                ggplot2::aes_string(label = "percent"),
                data     = data[data$prop > 5, ],
                position = ggplot2::position_dodge(width = 0.9),
                vjust    = 1.6,
                size     = 3.5,
                colour   = "#FFFFFF") +
            ggplot2::geom_text(
                ggplot2::aes_string(label = "percent", y = 0),
                data     = data[data$prop <= 5, ],
                position = ggplot2::position_dodge(width = 0.9),
                vjust    = 1.5,
                size     = 3.5,
                colour   = "#4D4D4D")
    }

    # Remove legend (if applicable)
    if (!legend) {
        gg <- gg + ggplot2::theme(legend.position = "none")
    }

    # Return graphics object
    return(gg)
}
