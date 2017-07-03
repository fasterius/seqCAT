#' Filter variants on specified criteria
#'
#' This is a function for filtering variants on various criteria. Variants with
#' a depth lower than 10 are removed, as are variants that don't have a
#' complete genotype for both samples. 
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

    # Find variants that have complete genotypes in both samples
    alleles = c('A1.input_1', 'A2.input_1', 
                'A1.input_2', 'A2.input_2')
        
    # Find rows with data for both alleles
    complete = row.names(
        subset(variants,rowSums(is.na(variants[, alleles])) == 0))

    # Remove non-complete genotypes
    variants = variants[complete, ]

    # Return the filtered variants
    return(variants)
}
