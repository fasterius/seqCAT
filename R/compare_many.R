#' @title Comparisons of many SNV profiles
#'
#' @description \link{compare_many} overlaps  and compares genotypes in many
#'  SNV profiles
#'
#' @details This is a function that compares all the combinations of the SNV
#' profiles input to it, either in a one-to-many or many-to-many manner. It
#' returns both a dataframe containing summary statistics for all unique
#' combinations and a list of dataframes with all the performed comparisons,
#' for easy re-use and downstream analyses of said comparisons.
#'
#' @export
#' @rdname compare_many
#' @param many List of SNV profiles to be compared
#' @param one Profile to be compared to all others
#' @param a Similarity score parameter a
#' @param b Similarity score parameter b
#' @return A list of summary statistics and comparisons
#'
#' @examples
#' # Load test data
#' data(test_profile_1)
#' data(test_profile_2)
#'
#' # Perform many-to-many comparisons
#' profiles <- list(test_profile_1, test_profile_2)
#' comparisons <- compare_many(profiles)
#'
#' # View aggregate similarities
#' \dontrun{comparisons[[1]])}
#'
#' # View data of first comparison
#' \dontrun{head(comparisons[[2]][[1]])}
compare_many <- function(many,
                         one  = NULL,
                         a    = 1,
                         b    = 5) {

    # Initialise objects to be returned
    similarities <- data.frame(sample_1         = character(),
                               sample_2         = character(),
                               overlaps         = numeric(),
                               matches          = numeric(),
                               concordance      = numeric(),
                               similarity_score = numeric(),
                               stringsAsFactors = FALSE)
    comparisons <- list()

    # One-to-many comparisons
    if (!is.null(one)) {

        # Find sample for one
        sample_one <- unique(one$sample)

        # Perform one's self-comparison
        comparison <- compare_profiles(one, one)
        similarities <- calculate_similarity(comparison,
                                             similarity = similarities,
                                             a          = a,
                                             b          = b)

        # Add comparison to collection
        comparisons[[length(comparisons) + 1]] <- comparison

        # Perform all one-to-many comparisons
        for (current in many) {


            # Find current's sample
            sample_current <- unique(current$sample)

            # Skip if current combination already exists in data
            exists <- similarities[similarities$sample_1 == sample_one &
                                   similarities$sample_2 == sample_current, ]
            if (nrow(exists) != 0) {
                next
            }

            # Compare and calculate
            comparison <- compare_profiles(one, current)
            similarities <- calculate_similarity(comparison,
                                                 similarity = similarities,
                                                 a          = a,
                                                 b          = b)

            # Add comparison to collection
            comparisons[[length(comparisons) + 1]] <- comparison
        }

    } else {

        # Many-to-many comparisons
        for (current_1 in many) {

            # Find current_1's sample
            sample_current_1 <- unique(current_1$sample)

            for (current_2 in many) {

                # Find current_2's sample
                sample_current_2 <- unique(current_2$sample)

                # Skip if current combination already exists in data
                exists <- similarities[
                    (similarities$sample_1 == sample_current_1 &
                     similarities$sample_2 == sample_current_2) |
                    (similarities$sample_1 == sample_current_2 &
                     similarities$sample_2 == sample_current_1), ]
                if (nrow(exists) != 0) {
                    next
                }

                # Compare and calculate
                comparison <- compare_profiles(current_1, current_2)
                similarities <- calculate_similarity(comparison,
                                                     similarity = similarities,
                                                     a          = a,
                                                     b          = b)

                # Add comparison to collection
                comparisons[[length(comparisons) + 1]] <- comparison
            }
        }
    }

    # Return list of results
    return(list(similarities, comparisons))
}
