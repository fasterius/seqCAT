#' Find matching variants in a data frame.
#'
#' This is a function for finding matching or mismatching variants across
#' two samples contained in a single data frame.
#'
#' A new column ("match") will be added to the data frame, which indicates if
#' each variant is a perfect match or not (i.e. "match" or "mismatch").
#'
#' @param data The dataframe to analyse.
#' @return A dataframe with an added "match" column
#' @examples
#' variant_matches(variants.df)

#' @export
variant_matches = function(data) {

    # Get sample names
    sample_1 = unique(data$sample.input_1)
    sample_2 = unique(data$sample.input_2)

    # Specify columns containing the alleles
    alleles = c('A1.input_1', 'A2.input_1', 'A1.input_2', 'A2.input_2')

    # Find variants with full genotypes in both samples
    idx.notna = row.names(subset(data, rowSums(is.na(data[, alleles])) == 0))

    # Check if there are overlapping variants
    if ( length(idx.notna) != 0 ) {  # At least one overlapping variant

        # Remove non-overlapping variants
        data = data[idx.notna, ]

        # Check for genotype matches
        data$match = 'mismatch'
        data.alleles = data[alleles]
        data.alleles$in.1 = paste(data.alleles[, 1], 
                                  data.alleles[, 2], sep=':')
        data.alleles$in.2 = paste(data.alleles[, 3], 
                                  data.alleles[, 4], sep=':')
        data.alleles$in.1.rev = paste(data.alleles[, 2], 
                                      data.alleles[, 1], sep=':')

        idx.match.1 = apply(data.alleles, 1, 
                            function(x) x['in.1'] %in% x['in.2'])
        idx.match.2 = apply(data.alleles, 1, 
                            function(x) x['in.1.rev'] %in% x['in.2'])
        data[idx.match.1, 'match'] = 'match'
        data[idx.match.2, 'match'] = 'match'

    } else {  # No overlapping variants

        # Re-introduce sample names
        data[1, 'sample.input_1'] = sample_1
        data[1, 'sample.input_2'] = sample_2

    }

    # Return the results
    return(data)

}
