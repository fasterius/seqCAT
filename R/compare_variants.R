#' Find matching variants in a data frame.
#'
#' This is a function for comparing matching or mismatching variants across
#' two samples contained in a single data frame.
#'
#' A new column ('match') will be added to the data frame, which indicates if
#' each variant's genotype is a match or not (i.e. 'match' or 'mismatch').
#'
#' @param overlaps The dataframe to analyse.
#' @return A dataframe with an added 'match' column
#' @examples
#' data(overlaps)
#' compare_variants(overlaps)

#' @export
compare_variants <- function(overlaps) {

    # Get sample names
    sample_1 <- unique(overlaps$sample_1)
    sample_2 <- unique(overlaps$sample_2)

    # Find overlapping variants with complete genotypes
    alleles <- paste(c("A1", "A1", "A2", "A2"),
                     c(sample_1, sample_2),
                     sep = ".")
    idx_notna <- row.names(
        overlaps[stats::complete.cases(overlaps[, alleles]), ])

    # Check for matches if there are overlapping variants
    if (length(idx_notna) != 0) {

        # Set all to "mismatch"
        overlaps$match <- "mismatch"

        # Construct alleles
        overlaps_alleles <- overlaps[alleles]
        overlaps_alleles$in_1 <- paste(overlaps_alleles[, 1],
                                   overlaps_alleles[, 3],
                                   sep = ":")
        overlaps_alleles$in_2 <- paste(overlaps_alleles[, 2],
                                   overlaps_alleles[, 4],
                                   sep = ":")
        overlaps_alleles$in_1_rev <- paste(overlaps_alleles[, 3],
                                       overlaps_alleles[, 1],
                                       sep = ":")

        # Check and set matching genotypes as appropriate
        idx_match_1 <- apply(overlaps_alleles, 1,
                            function(x) x["in_1"] %in% x["in_2"])
        idx_match_2 <- apply(overlaps_alleles, 1,
                            function(x) x["in_1_rev"] %in% x["in_2"])
        overlaps[idx_match_1, "match"] <- "match"
        overlaps[idx_match_2, "match"] <- "match"

    } else {

        # Add empty match column if no overlapping variants are found
        overlaps$match <- NA
    }

    # Return the results
    return(overlaps)

}
