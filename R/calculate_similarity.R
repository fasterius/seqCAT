#' @title SNV profile similarity calculations
#'
#' @description \link{calculate_similarity} returns the summary statistics for
#'  SNV profile comparisons
#'
#' @details
#' This function calculates various summary statistics and sample similarities
#' for a given profile comparison dataframe. It returns a small dataframe with
#' the overall similarity score (whose parameters `a` and `b` can be adjusted
#' in the function call), total SNV overlaps, the concordance of the overlaps
#' and the sample names in question. This dataframe can also be given to the 
#' function, in which case it will simply add another row for the current
#' samples, facilitating downstream aggregate analyses.
#'
#' @export
#' @rdname calculate_similarity
#' @param overlaps The input SNV overlaps dataframe
#' @param a Similarity score parameter a
#' @param b Similarity score parameter b
#' @param similarity Optional dataframe to add results to
#' @return A dataframe with summary statistics
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
calculate_similarity <- function(overlaps,
                                 similarity = NULL,
                                 a          = 1,
                                 b          = 5) {

    # Check similarity dataframe structure (if applicable)
    if (!is.null(similarity)) {

        # Check that similarity is a dataframe
        if (class(similarity) != "data.frame") {

            stop("supplied similarity object is not a dataframe")

        }

        # Correct dataframe structure
        correct_names <- c("sample_1",
                           "sample_2",
                           "overlaps",
                           "matches",
                           "concordance",
                           "similarity_score")

        if (!identical(names(similarity), correct_names)) {

            stop(paste0("supplied similarity dataframe does not have the ",
                        "correct structure"))

        }
    }

    # Get samples
    sample_1 <- unique(overlaps$sample_1)
    sample_2 <- unique(overlaps$sample_2)

    # Get number of overlaps and matches
    n_total <- nrow(overlaps)
    n_matches <- nrow(overlaps[overlaps$match == "match", ])

    # Calculate concordance
    concordance <- round(n_matches / n_total * 100, 1)

    # Calculate similarity score
    score <- round( (n_matches + a) / (n_total + a + b) * 100, 1)

    # Create results dataframe
    results <- data.frame(sample_1         = sample_1,
                          sample_2         = sample_2,
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
