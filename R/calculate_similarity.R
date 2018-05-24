#' @title SNV profile similarity calculations
#'
#' @description Calculate the similarity statistics for SNV profile
#'  comparisons.
#'
#' @details
#' This function calculates various summary statistics and sample similarities
#' for a given profile comparison dataframe. It returns a small dataframe with
#' the overall similarity score (whose parameters `a` and `b` can be adjusted
#' in the function call), total SNV data, the concordance of the data
#' and the sample names in question. This dataframe can also be given to the 
#' function, in which case it will simply add another row for the current
#' samples, facilitating downstream aggregate analyses.
#'
#' @export
#' @rdname calculate_similarity
#' @param data The input SNV data dataframe.
#' @param a Similarity score parameter a (integer).
#' @param b Similarity score parameter b (integer).
#' @param similarity Optional dataframe to add results to.
#' @return A dataframe with summary statistics.
#'
#' @examples
#' # Load test data
#' data(test_comparison)
#'
#' # Calculate similarities
#' similarity <- calculate_similarity(test_comparison)
#'
#' # Add another row of summary statistics
#' calculate_similarity(test_comparison, similarity = similarity)
calculate_similarity <- function(data,
                                 similarity = NULL,
                                 a          = 1,
                                 b          = 5) {

    # Check similarity dataframe structure (if applicable)
    if (!is.null(similarity)) {

        # Check that similarity is a dataframe
        if (!methods::is(similarity, "data.frame")) {
            stop("supplied similarity object is not a dataframe")
        }

        # Correct dataframe structure
        correct_names <- c("sample_1",
                           "sample_2",
                           "variants_1",
                           "variants_2",
                           "overlaps",
                           "matches",
                           "concordance",
                           "similarity_score")

        if (!identical(names(similarity), correct_names)) {
            stop("supplied similarity dataframe does not have the ",
                 "correct structure")
        }
    }

    # Get samples
    sample_1 <- unique(data$sample_1)
    sample_2 <- unique(data$sample_2)

    # Get number of overlaps and matches
    n_total <- nrow(data[data$match == "match" |
                         data$match == "mismatch", ])
    n_matches <- nrow(data[data$match == "match", ])

    # Get number of non-overlaps, if present
    n_sample_1 <- nrow(data) -
        nrow(data[data$match == paste0(sample_2, "_only"), ])
    n_sample_2 <- nrow(data) -
        nrow(data[data$match == paste0(sample_1, "_only"), ])

    # Calculate concordance
    concordance <- round(n_matches / n_total * 100, 1)

    # Calculate similarity score
    score <- round( (n_matches + a) / (n_total + a + b) * 100, 1)

    # Create results dataframe
    results <- data.frame(sample_1         = sample_1,
                          sample_2         = sample_2,
                          variants_1       = n_sample_1,
                          variants_2       = n_sample_2,
                          overlaps         = n_total,
                          matches          = n_matches,
                          concordance      = concordance,
                          similarity_score = score,
                          stringsAsFactors = FALSE)

    # Combine with existing results (if applicable)
    if (!is.null(similarity)) {
        results <- rbind(similarity, results)
    }

    # Return results dataframe
    return(results)
}
