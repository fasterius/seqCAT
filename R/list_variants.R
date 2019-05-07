#' @title List known variants
#'
#' @description List known variants present in SNV profiles
#'
#' @details This is a function for listing known variants present in SNV
#'  profiles. Input is a list of profiles and a dataframe of known variants,
#'  containing at least the genomic locations ("chr" and "pos"). Any additional
#'  columns will be retained.
#'
#' @export
#' @param profiles The SNV profiles to analyse (list)
#' @param known_variants The known variants to look for (dataframe)
#' @return A dataframe containing the known variant genotypes in each profile.
#'
#' @examples
#' # Load test data
#' data(test_profile_1)
#' data(test_profile_2)
#'
#' # Create some variants to analyse
#' known_variants <- data.frame(chr = 1, pos = 16229, gene = "DDX11L1")
#'
#' # List the known variants in each profile
#' profiles <- list(test_profile_1, test_profile_2)
#' known_variants <- list_variants(profiles, known_variants)
list_variants <- function(profiles, known_variants) {

    # Test that `known_variants` is a dataframe
    if (!methods::is(known_variants, "data.frame")) {
        stop("`known_variants` input is a not a dataframe.")
    }

    # Test that known_variants contains genomic locations
    if (!("chr" %in% names(known_variants))) {
        stop("`known_variants` dataframe does not contain 'chr' column.")
    } else if (!("pos" %in% names(known_variants))) {
        stop("`known_variants` dataframe does not contain 'pos' column.")
    }

    # Create GRanges object from known variants
    known_gr <- GenomicRanges::makeGRangesFromDataFrame(known_variants,
        keep.extra.columns      = TRUE,
        ignore.strand           = TRUE,
        seqinfo                 = NULL,
        seqnames.field          = "chr",
        start.field             = "pos",
        end.field               = "pos",
        starts.in.df.are.0based = FALSE)

    # Loop through each profile and compare to known variants
    for (profile in profiles) {

        # Convert profile to GRanges object
        profile <- convert_to_gr(profile)

        # Compare to known variants
        sample <- unique(profile$sample)
        comp <- S4Vectors::intersect(profile, known_gr)
        comp <- add_metadata(comp, profile, "")

        # Add zeroes if no known variants are present
        if (length(comp) == 0) {
            known_variants[[sample]] <- 0
            next
        }

        # Convert to dataframe
        comp <- GenomicRanges::as.data.frame(comp)

        # Remove duplicated positions
        comp <- comp[!duplicated(comp[, c('seqnames', 'start')]), ]

        # Combine alleles
        comp[[sample]] <- paste0(comp$A1, "/", comp$A2)

        # Merge with known variants dataframe
        results <- comp[c("seqnames", "start", sample)]
        known_variants <- merge(known_variants,
                                results,
                                all.x           = TRUE,
                                by.x            = c("chr", "pos"),
                                by.y            = c("seqnames", "start"))
    }

    # Replace NAs with zeroes
    known_variants[is.na(known_variants)] <- 0

    # Return collated results
    return(known_variants)
}
