#' Plot SNV impact distributions for a binary comparison
#'
#' This function creates publication-ready plots of the impact distribution
#' from a binary dataset comparison across the matched/mismatched SNVs.
#'
#' @export
#' @rdname plot_impacts
#' @param comparison The profile comparison to be plotted
#' @param annotate Annotate each category
#' @param legend Show/hide legend
#' @param palette Colour palette for filling of bars
#' @examples
#' data(test_comparison)
#' plot_impacts(test_comparison)
plot_impacts <- function(comparison,
                         annotate = TRUE,
                         legend = TRUE,
                         palette = c("#0D2D59", "#1954A6")) {

    # Factorise match and impact columns
    comparison$match <- factor(comparison$match,
                               levels = c("match", "mismatch"))
    comparison$impact <-
        factor(comparison$impact,
               levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))

    # Calculate impact distribution
    # impact <- comparison %>%
    impact <- dplyr::group_by(comparison, match, impact)
    impact <- dplyr::summarise(impact, count = n())
    impact <- dplyr::mutate(impact, prop = count / sum(count) * 100)

    # Create plot object
    gg <- ggplot2::ggplot(impact, ggplot2::aes(x    = impact,
                                               y    = prop,
                                               fill = match)) +
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
        impact$percent <- paste0(format(impact$prop, nsmall = 1,
                                        digits = 1), " %")

        # Add annotations to plot
        gg <- gg +
            ggplot2::geom_text(
                ggplot2::aes(label = count),
                position           = ggplot2::position_dodge(width = 0.9),
                vjust              = -0.5,
                size               = 2.5,
                colour             = "#000000") +
            ggplot2::geom_text(
                data               = impact[impact$prop > 5, ],
                ggplot2::aes(label = percent),
                position           = ggplot2::position_dodge(width = 0.9),
                vjust              = 1.6,
                size               = 3.5,
                colour             = "#FFFFFF") +
            ggplot2::geom_text(
                data               = impact[impact$prop <= 5, ],
                ggplot2::aes(label = percent),
                position           = ggplot2::position_dodge(width = 0.9),
                vjust              = 1.6,
                size               = 3.5,
                colour             = "#4D4D4D")
    }

    # Remove legend (if applicable)
    if (!legend) {
        gg <- gg + ggplot2::theme(legend.position = "none")
    }

    # Return graphics object
    return(gg)
}
