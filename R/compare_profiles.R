#' @title Binary SNV profile comparisons
#'
#' @description Overlap and compare genotypes in two SNV profiles.
#'
#' @details This is a function for finding overlapping variants in two
#' different SNV profiles (stored as GenomicRanges objects), followed by
#' comparing the genotypes of the overlapping variants. The "compare_overlaps"
#' function calls the "add_metadata" function twice in succession in order to
#' merge the metadata for the two profiles (supplied as GRanges objects),
#' returns the results as a dataframe, compares the genotypes of the
#' overlapping variants using the "compare_genotypes" function and, finally,
#' returns the final dataframe with all variant overlaps and their similarity.

#' @export
#' @rdname compare_profiles
#' @param profile_1 The first SNV profile (GRanges object).
#' @param profile_2 The second SNV profile (GRanges object).
#' @param mode Merge profiles using "union" or "intersection" (character).
#' @return A dataframe.
#' 
#' @examples
#' # Load test data
#' data(test_profile_1)
#' data(test_profile_2)
#'
#' # Compare the two profiles
#' comparison <- compare_profiles(test_profile_1, test_profile_2)
compare_profiles <- function(profile_1,
                             profile_2,
                             mode = "intersection") {

    # Find sample names
    sample_1 <- unique(profile_1$sample)
    sample_2 <- unique(profile_2$sample)

    # Message
    message("Comparing ", sample_1, " and ", sample_2, " ...")

    # Find the overlaps of all ranges in both objects
    if (tolower(mode) == "union") {
        data_gr <- S4Vectors::union(profile_1, profile_2)
    } else if (tolower(mode) == "intersection") {
        data_gr <- S4Vectors::intersect(profile_1, profile_2)
    } else {
        stop("`mode` must be either 'union' or 'intersection'")
    }

    # Add metadata from both objects to the union object
    data_gr <- add_metadata(data_gr,
                            profile_1,
                            paste0(".", sample_1))
    data_gr <- add_metadata(data_gr,
                            profile_2,
                            paste0(".", sample_2))
    # Convert to data frame
    data <- GenomicRanges::as.data.frame(data_gr)

    # Compare genotypes in each overlapping SNV
    data <- compare_genotypes(data, sample_1, sample_2)

    # Collate metadata columns
    data <- collate_metadata(data, sample_1, sample_2)

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
compare_genotypes <- function(data, sample_1, sample_2) {

    # Remove non-complete variants
    alleles <- paste(c("A1", "A1", "A2", "A2"),
                     c(sample_1, sample_2),
                     sep = ".")
    data <- data[rowSums(is.na(data[, alleles])) == 0 |
                 rowSums(is.na(data[, alleles])) == 2, ]

    # Check if there are no variants in the data
    if (nrow(data) == 0) {

        # Manually add sample names and status
        data[1, "sample_1"] <- sample_1
        data[1, "sample_2"] <- sample_2
        data$match <- "no overlaps"

        # Return empty data
        return(data)
    } 

    # Add sample names
    data$sample_1 <- sample_1
    data$sample_2 <- sample_2

    # Set all to "mismatch"
    data$match <- "mismatch"

    # Construct alleles
    data_alleles <- data[alleles]
    data_alleles$in_1 <- paste(data_alleles[, 1], data_alleles[, 3], sep = ":")
    data_alleles$in_2 <- paste(data_alleles[, 2], data_alleles[, 4], sep = ":")
    data_alleles$rev <- paste(data_alleles[, 3], data_alleles[, 1], sep = ":")

    # Check and set matching genotypes as appropriate
    idx_match_1 <- apply(data_alleles, 1, function(x) x["in_1"] %in% x["in_2"])
    idx_match_2 <- apply(data_alleles, 1, function(x) x["rev"] %in% x["in_2"])
    data[idx_match_1, "match"] <- "match"
    data[idx_match_2, "match"] <- "match"

    # Check and set non-overlapping variants as appropriate
    alleles_1 <- paste(c("A1", "A2"), sample_1, sep = ".")
    alleles_2 <- paste(c("A1", "A2"), sample_2, sep = ".")
    data[rowSums(is.na(data[, alleles_1])) == 0 &
         rowSums(is.na(data[, alleles_2])) != 0, "match"] <-
             paste0(sample_1, "_only")
    data[rowSums(is.na(data[, alleles_2])) == 0 &
         rowSums(is.na(data[, alleles_1])) != 0, "match"] <-
             paste0(sample_2, "_only")

    # Return the results
    return(data)
}

# Function for collating metadata columns from both samples
collate_metadata <- function(data, sample_1, sample_2) {

    # Remove redundant sample name columns
    data <- dplyr::select_(data,
                           paste0("-sample.", sample_1),
                           paste0("-sample.", sample_2))

    # Find common, sample-specific metadata columns
    s1 <- paste0("\\.", sample_1)
    s2 <- paste0("\\.", sample_2)
    mcols_sample_1 <- gsub(s1, "", grep(s1, names(data), value = TRUE))
    mcols_sample_2 <- gsub(s2, "", grep(s2, names(data), value = TRUE))
    mcols <- mcols_sample_1[mcols_sample_1 %in% mcols_sample_2]

    # Remove data-specific metadata columns
    mcols <- grep("DP|AD1|AD2|A1|A2|warnings",
                  mcols,
                  value = TRUE, invert = TRUE)

    # Loop over metadata columns and merge as applicable
    for (mcol in mcols) {

        # Current metadata column name
        mcol_s1 <- paste0(mcol, ".", sample_1)
        mcol_s2 <- paste0(mcol, ".", sample_2)

        # Create merged metadata column as appropriate
        if (!(all(is.na(data[[mcol_s1]]), is.na(data[[mcol_s2]])))) {

            # Get index for identical metadata
            idx <- data[[mcol_s1]] == data[[mcol_s2]]
            idx_na <- is.na(idx)
            idx[is.na(idx)] <- TRUE

            # Create merged metadata
            data[idx & !is.na(idx), mcol] <- data[idx, mcol_s1]
            data[is.na(data[idx, mcol]), mcol] <-
                data[is.na(data[idx, mcol]), mcol_s2]
            data[idx & is.na(idx), mcol] <- paste0(data[is.na(idx), mcol_s1],
                                                   data[is.na(idx), mcol_s2])
            data[!idx & !is.na(!idx), mcol] <-
                paste0("[", data[!idx, mcol_s1],
                       ",", data[!idx, mcol_s2], "]")
        }

        # Remove old metadata columns
        data <- dplyr::select_(data,
                               paste0("-", mcol_s1),
                               paste0("-", mcol_s2))
    }

    # Delete unneccesary columns
    data <- dplyr::select_(data,
                          "-end",
                          "-width",
                          "-strand")

    # Rename columns
    names(data)[c(1, 2)] <- c("chr", "pos")

    # Fix column names for COSMIC comparisons
    cols <- c("rsID", "ENSGID", "ENSTID", "impact",
              "effect", "feature", "biotype")
    for (col in cols) {
        names(data)[grep(col, names(data))] <- col
    }

    # Re-order columns
    with_sample <- grep("\\.", names(data), value = TRUE)
    without_sample <- grep("\\.", names(data), value = TRUE, invert = TRUE)
    firsts <- c("chr", "pos", "sample_1", "sample_2", "match")
    without_sample <- setdiff(without_sample, firsts)
    data <- data[c(firsts, without_sample, with_sample)]

    # Remove redundant columns for comparisons without overlaps
    if (nrow(data) == 1 & data[1, "match"] == "no overlaps") {
        data <- data[c("sample_1", "sample_2", "match")]
    }

    # Remove "<NA>" strings
    data[is.na(data)] <- ""

    # Return final data
    return(data)
}
