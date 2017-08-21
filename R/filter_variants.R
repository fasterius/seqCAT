#' Filter variants on specified criteria
#'
#' This is a function for filtering variants on sequencing depth. Variants with
#' a depth lower than 10 are removed by default, but can be changed in the
#' function call.
#'
#' @param overlaps The data frame containing the variant data to be filtered.
#' @param filter_depth Threshold for variant depth (default 10)
#' @return A data frame containing the filtered variants.
#' @examples
#' data(test_overlaps)
#' filter_variants(test_overlaps)
#' filter_variants(test_overlaps, filter_depth = 20)

#' @export
filter_variants <- function(overlaps, filter_depth = 10) {

    # Find sample names
    sample_1 <- unique(overlaps$sample_1)
    sample_2 <- unique(overlaps$sample_2)

    # Filter on sequencing depth
    overlaps <- overlaps[overlaps[[paste0("DP.", sample_1)]] >= filter_depth &
                         overlaps[[paste0("DP.", sample_2)]] >= filter_depth, ]

    # Return the filtered overlaps
    return(overlaps)
}
