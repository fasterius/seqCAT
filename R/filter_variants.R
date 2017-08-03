#' Filter variants on specified criteria
#'
#' This is a function for filtering variants on sequencing depth. Variants with
#' a depth lower than 10 are removed by default, but can be changed in the
#' function call.
#'
#' @param variants The data frame containing the variant data to be filtered.
#' @param filter_depth Threshold for variant depth (default 10)
#' @return A data frame containing the filtered variants.
#' @examples
#' filter_variants(data)
#' filter_variants(data, filter_depth=20)

#' @export
filter_variants = function(variants, filter_depth=10) {

    # Filter on sequencing depth
    variants = variants[variants$DP.input_1 >= filter_depth & 
                        variants$DP.input_2 >= filter_depth, ]

    # Return the filtered variants
    return(variants)
}
