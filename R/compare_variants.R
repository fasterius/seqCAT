#' Find matching variants in a data frame.
#'
#' This is a function for comparing matching or mismatching variants across
#' two samples contained in a single data frame.
#'
#' A new column ("match") will be added to the data frame, which indicates if
#' each variant's genotype is a match or not (i.e. "match" or "mismatch").
#'
#' @param data The dataframe to analyse.
#' @return A dataframe with an added "match" column
#' @examples
#' compare_variants(variants.df)

#' @export
compare_variants = function(data) {

    # Get sample names
    sample_1 = unique(data$sample_1)
    sample_2 = unique(data$sample_2)

    # Find for overlapping variants with complete genotypes
    alleles = paste(c("A1", "A2", "A1", "A2"), c(sample_1, sample_2), sep=".")
    idx.notna = row.names(data[complete.cases(data[, alleles]), ])

    # Check for matches if there are any overlapping variants
    if ( length(idx.notna) != 0 ) {

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
    }

    # Return the results
    return(data)

}
