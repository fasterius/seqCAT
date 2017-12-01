#' @title Plot known variants list
#'
#' @description Plot a genotype grid from a list of known variants
#'
#' @details This function creates publication-ready plots from lists of known
#'  variants, taking a dataframe containing all the genotypes (on "A1/A2"
#'  format) for each sample (columns) and variant (row names).
#'
#' @export
#' @rdname plot_variant_list
#' @param variant_list The data containing the variants (dataframe)
#' @param legend Show a legend for the genotype colours (boolean)
#' @param palette Nucleotide colour palette (4-element character vector)
#' @return A ggplot2 graphical object.
#'
#' @examples
#' # Load test variant list
#' data(test_variant_list)
#'
#' # Plot each variant's genotype per sample
#' genotype_grid <- plot_variant_list(test_variant_list)
plot_variant_list <- function(variant_list, 
                              legend        = TRUE,
                              palette       = c("#4e8ce4", "#a6c6f2",
                                                "#999999", "#cccccc")) {
    
    # Convert to wide format
    variant_list$variant <- row.names(variant_list)
    variants <- tidyr::gather_(data       = variant_list,
                               key        = "sample",
                               value      = "genotype",
                               setdiff(names(variant_list), "variant"))
    
    # add ones to each non-zero cell
    variants[variants$genotype == 0, "genotype"] <- NA
    
    # Separate genotypes into alleles
    variants <- tidyr::separate_(data   = variants,
                                 col    = "genotype",
                                 sep    = "/",
                                 into   = c("A1", "A2"),
                                 remove = TRUE)

    # Factorise alleles
    allele_levels <- c("A", "T", "C", "G")
    variants$A1 <- factor(variants$A1, levels = allele_levels)
    variants$A2 <- factor(variants$A2, levels = allele_levels)
    
    # Plot
    genotype_grid <- ggplot2::ggplot(variants,
                                     ggplot2::aes_string(x      = "variant",
                                                         y      = "sample",
                                                         fill   = "A1",
                                                         colour = "A2")) +
        ggplot2::geom_tile(size = 0.3) +
        ggplot2::coord_equal() +
        ggplot2::theme(axis.ticks       = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       axis.text.x      = ggplot2::element_text(angle = 90,
                                                                hjust = 1,
                                                                vjust = 0.5)) +
        ggplot2::labs(x    = NULL,
                      y    = NULL,
                      fill = NULL) +
        ggplot2::scale_fill_manual(breaks = c("A", "T", "C", "G"),
                                   values = palette, 
                                   drop = FALSE) +
        ggplot2::scale_colour_manual(breaks = c("A", "T", "C", "G"),
                                     values = palette,
                                     drop = FALSE) +
        ggplot2::guides(colour = FALSE)

    # Add polygons for A2
    genotype_grid <- add_polygons(genotype_grid, allele_levels, palette)

    # Add legend (if applicable)
    if (!legend) {
        genotype_grid <- genotype_grid +
            ggplot2::theme(legend.position = "none")
    }

    # Return plot
    return(genotype_grid)
}

# Function for adding allele colours
add_polygons <- function(genotype_grid, allele_levels, palette) {

    # Build graphical data from existing genotype_grid
    gg_object <- ggplot2::ggplot_build(genotype_grid)
    gg_data <- gg_object$data[[1]]

    # Reset tile edge colour to black
    genotype_grid <- genotype_grid + ggplot2::geom_tile(colour = "black")
    
    # Subset for coloured tiles
    gg_data <- gg_data[!is.na(gg_data$fill), ]

    # Get groups per tile
    groups <- factor(gg_data$group)

    # Get colours
    colours <- data.frame(group  = groups,
                          colour = gg_data$colour)
    colours$colour <- as.character(colours$colour)

    # Reverse-map colours to alleles
    colours[colours$colour == palette[1], "A1"] <- "A"
    colours[colours$colour == palette[2], "A1"] <- "T"
    colours[colours$colour == palette[3], "A1"] <- "C"
    colours[colours$colour == palette[4], "A1"] <- "G"
    colours$A1 <- factor(colours$A1, levels = allele_levels)
    colours$colour <- NULL

    # Get positions
    x <- data.frame(x1 = gg_data$xmin,
                    x2 = gg_data$xmin,
                    x3 = gg_data$xmax,
                    x4 = gg_data$xmin)
    y <- data.frame(y1 = gg_data$ymin,
                    y2 = gg_data$ymax,
                    y3 = gg_data$ymax,
                    y4 = gg_data$ymin)
    positions <- data.frame(group = rep(groups, each = 4),
                            x     = as.vector(t(x)),
                            y     = as.vector(t(y)))

    # Merge groups/colours and groups/positions
    polygon_data <- merge(colours, positions, by = "group")

    # Add polygons to genotype_grid
    genotype_grid <- genotype_grid +
        ggplot2::geom_polygon(data   = polygon_data,
                              colour = "black",
                              size   = 0.3,
                              ggplot2::aes_string(x      = "x",
                                                  y      = "y",
                                                  fill   = "A1",
                                                  group  = "group"))

    # Return genotype_grid with polygons
    return(genotype_grid)
}
