#' Filter variants on specified criteria
#'
#' This is a function for filtering variants on various criteria. Variants with
#' a depth lower than 10 are removed, as are variants that don't have a
#' complete genotype for both samples (by default).
#'
#' @param variants The data frame containing the variant data to be filtered.
#' @param overlaps Boolean that governs filtering on complete genotypes
#'     (enabled by default)
#' @return A data frame containing the filtered variants.
#' @examples
#' filter_variants(data)
#' filter_variants(data, overlaps = FALSE)

filter_variants = function(variants, overlaps = TRUE) {

    # Filter on sequencing depth for each variant set
    variants = variants[variants$DP.input_1 >= 10 & 
                        variants$DP.input_2 >= 10, ]
    
    # Filter on the "filter" column from GATK
    variants = variants[variants$filter.input_1 == 'None' &
                        variants$filter.input_2 == 'None', ]

    # Filter to only include overlaps (default)
    if ( overlaps ) {
        
        # Find variants that have complete genotypes in both samples
        alleles = c('allele_1.input_1', 'allele_2.input_1', 
                    'allele_1.input_2', 'allele_2.input_2')
        
        # Find rows with data for both alleles
        complete = row.names(subset(variants, 
                                    rowSums(is.na(variants[, alleles])) == 0))

        # Remove non-complete genotypes
        variants = variants[complete, ]
    }

    # Return the filtered variants
    return(variants)
}
