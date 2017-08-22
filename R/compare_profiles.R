#' Overlap and compare genotypes in two SNV profiles.
#'
#' This is a function for finding overlapping variants in two different SNV
#' profiles (stored as GenomicRanges objects), followed by comparing the
#' genotypes of the overlapping variants. The "compare_overlaps" function calls
#' the "add_metadata" function twice in succession in order to merge the
#' metadata for the two profiles (supplied as GRanges objects), returns the
#' results as a dataframe, compares the genotypes of the overlapping variants
#' using the "compare_genotypes" function and, finally, returns the final
#' dataframe with all variant overlaps and their similarity.

#' @export
#' @rdname compare_profiles
#' @param profile_1 The first variant GRanges object.
#' @param profile_2 The second variant GRanges object.
#' @param sample_1 Name of the first sample.
#' @param sample_2 Name of the second sample.
#' @return A data frame.
#' @examples
#' data(test_profile_1)
#' data(test_profile_2)
#' compare_profiles(test_profile_1, test_profile_2, "sample1", "sample2")
compare_profiles <- function(profile_1,
                             profile_2,
                             sample_1 = "sample_1",
                             sample_2 = "sample_2") {

    # Find the intersection of all ranges in both objects
    intersect_gr <- S4Vectors::intersect(profile_1, profile_2)

    # Add metadata from both objects to the union object
    intersect_gr <- add_metadata(intersect_gr,
                                 profile_1,
                                 paste0(".", sample_1))
    intersect_gr <- add_metadata(intersect_gr,
                                 profile_2,
                                 paste0(".", sample_2))

    # Convert to data frame
    data <- GenomicRanges::as.data.frame(intersect_gr)

    # Remove non-complete variants
    alleles <- paste(c("A1", "A1", "A2", "A2"),
                     c(sample_1, sample_2),
                     sep = ".")
    data <- data[stats::complete.cases(data[, alleles]), ]

    # Add empty data frame with sample names if no variants overlap
    if (nrow(data) == 0) {
        data[1, "sample_1"] <- sample_1
        data[1, "sample_2"] <- sample_2
    } else {
        data$sample_1 <- sample_1
        data$sample_2 <- sample_2
    }

    # Remove redundant sample name columns
    redundant <- paste0("sample.", c(sample_1, sample_2))
    data <- data[!names(data) %in% redundant]

    # Compare genotypes in each overlapping SNV
    data <- compare_genotypes(data)

    # Collate metadata columns
    data <- collate_metadata(data)

    # Return the final data frame
    return(data)
}

# Function for adding metadata from a GRanges <subject> to a GRanges <query>,
# while adding <column_suffix> to the added metadata columns.
add_metadata <- function(query,
                         subject,
                         column_suffix) {

    # Find overlapping ranges
    hits <- IRanges::findOverlaps(query, subject)

    for (column in names(S4Vectors::mcols(subject))) {

        # Create empty metadata column to be filled
        S4Vectors::mcols(query)[paste(column, column_suffix, sep = "")] <- NA

        # Convert DNAStringSet / DNAStringSetList columns to character vectors
        if (class(S4Vectors::mcols(subject)[[column]])[1] == "DNAStringSet") {
          S4Vectors::mcols(subject)[column] <-
              as.character(S4Vectors::mcols(subject)[[column]])
        } else if (class(S4Vectors::mcols(subject)[[column]])[1] ==
                   "DNAStringSetList") {
          S4Vectors::mcols(subject)[column] <-
            S4Vectors::unstrsplit(IRanges::CharacterList(
                S4Vectors::mcols(subject)[[column]]))
        }

    # Add subject metadata to query
    S4Vectors::mcols(query)[S4Vectors::queryHits(hits),
                            paste(column, column_suffix, sep = "")] <-
      S4Vectors::mcols(subject)[S4Vectors::subjectHits(hits), column]
    }
    return(query)
}

# Function for comparing genotypes in overlapping SNVs
compare_genotypes <- function(overlaps) {

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

# Function for collating metadata columns from both samples
collate_metadata <- function(data) {

    # Find sample-specific metadata columns
    mcols <- grep("\\.sample_1", names(data),
                  value = TRUE)
    mcols <- grep("DP|AD1|AD2|A1|A2|warnings", mcols,
                  value = TRUE, invert = TRUE)

    # Loop over metadata columns and merge as applicable
    for (mcol_1 in mcols) {

        # Current metadata column name
        mcol <- gsub("\\.sample_1", "", mcol_1)
        mcol_2 <- gsub("\\.sample_1", "\\.sample_2", mcol_1)

        # Create merged metadata column as appropriate
        if (!(all(is.na(data[[mcol_1]]), is.na(data[[mcol_2]])))) {

            data[data[[mcol_1]] == data[[mcol_2]], mcol] <- data[[mcol_1]]
            data[data[[mcol_1]] != data[[mcol_2]], mcol] <-
                paste0("[", data[[mcol_1]], ",", data[[mcol_2]], "]")
        }

        # Remove old metadata columns
        data <- dplyr::select_(data,
                               paste0("-", mcol_1),
                               paste0("-", mcol_2))
    }

    # Delete unneccesary columns
    data <- dplyr::select_(data,
                          "-end",
                          "-width",
                          "-strand")

    # Rename columns
    names(data)[1:2] <- c("chr", "pos")

    # Return final data
    return(data)
}
