#' Find matching variants in a data frame.
#'
#' This is a function for comparing matching or mismatching variants across
#' two samples contained in a single data frame.
#'
#' A new column ('match') will be added to the data frame, which indicates if
#' each variant's genotype is a match or not (i.e. 'match' or 'mismatch').
#'
#' @param data The dataframe to analyse.
#' @return A dataframe with an added 'match' column
#' @examples
#' compare_variants(variants.df)

#' @export
compare_variants <- function(data) {

    # Get sample names
    sample_1 <- unique(data$sample_1)
    sample_2 <- unique(data$sample_2)

    # Find overlapping variants with complete genotypes
    alleles <- paste(c("A1", "A1", "A2", "A2"),
                     c(sample_1, sample_2),
                     sep = ".")
    idx_notna <- row.names(data[complete.cases(data[, alleles]), ])

    # Check for matches if there are overlapping variants
    if (length(idx_notna) != 0) {

        # Set all to "mismatch"
        data$match <- "mismatch"

        # Construct alleles
        data_alleles <- data[alleles]
        data_alleles$in_1 <- paste(data_alleles[, 1],
                                   data_alleles[, 3],
                                   sep = ":")
        data_alleles$in_2 <- paste(data_alleles[, 2],
                                   data_alleles[, 4],
                                   sep = ":")
        data_alleles$in_1_rev <- paste(data_alleles[, 3],
                                       data_alleles[, 1],
                                       sep = ":")

        # Check and set matching genotypes as appropriate
        idx_match_1 <- apply(data_alleles, 1,
                            function(x) x["in_1"] %in% x["in_2"])
        idx_match_2 <- apply(data_alleles, 1,
                            function(x) x["in_1_rev"] %in% x["in_2"])
        data[idx_match_1, "match"] <- "match"
        data[idx_match_2, "match"] <- "match"

    } else {

        # Add empty match column if no overlapping variants are found
        data$match <- NA
    }

    # Return the results
    return(data)

}
