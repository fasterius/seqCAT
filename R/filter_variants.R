#' @title Variant filtering
#'
#' @description Filter variants on sequencing depth.
#'
#' @details This is a function for filtering variants on sequencing depth.
#' Variants with a depth lower than 10 are removed by default, but can be
#' changed in the function call.
#'
#' @export
#' @param data The dataframe containing the variant data to be filtered.
#' @param min_depth Threshold for variant depth (integer; default 10).
#' @return A data frame containing the filtered variants.
#'
#' @examples
#' # Load test comparisons
#' data(test_comparison)
#'
#' # Filter variants
#' filt_1 <- filter_variants(test_comparison)
#' filt_2 <- filter_variants(test_comparison, min_depth = 20)
filter_variants <- function(data,
                            min_depth = 10) {

    # Find sample names
    sample_1 <- unique(data$sample_1)
    sample_2 <- unique(data$sample_2)

    # Filter on sequencing depth
    col_1 <- paste0("DP.", sample_1)
    col_2 <- paste0("DP.", sample_2)
    data <- data[(data[[col_1]] >= min_depth | is.null(data[[col_1]])) &
                 (data[[col_2]] >= min_depth | is.null(data[[col_2]])), ]

    # Return the filtered data
    return(data)
}
