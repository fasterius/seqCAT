#' Find overlapping variants in two GenomicRanges objects.
#'
#' This are functions for finding overlapping variants in two different
#' GenomicRanges objects, returning another GRanges object.
#'
#' The add_metadata function is a function for adding metadata (i.e. any
#' column that is not the "seqnames", "start" or "end" fields in a
#' GenomicRanges object) from the "subject" GRanges object to the "query" 
#' GRanges object. The overlap_variants function, on the other hand, is a
#' wrapper function that calls add_metadata twice in succession, to add named
#' metadata to the intersection of variants between the two datasets.
#'
#' @param object_1 The first variant GRanges object.
#' @param object_2 The second variant GRanges object.
#' @param sample_1 Name of the first sample.
#' @param sample_2 Name of the second sample.
#' @return A GRanges object.
#' @examples
#' overlap_variants(data_first, data_second, "first_sample", "second_sample")

#' @rdname overlap_variants
add_metadata = function(query,
                        subject,
                        column_suffix) {

    # Find overlapping ranges
    hits = IRanges::findOverlaps(query, subject)

    for ( column in names(S4Vectors::mcols(subject)) ) {

        # Create empty metadata column to be filled
        S4Vectors::mcols(query)[paste(column, column_suffix, sep='')] = NA

        # Convert DNAStringSet / DNAStringSetList columns to character vectors
        if (class(S4Vectors::mcols(subject)[[column]])[1] == 'DNAStringSet') {
          S4Vectors::mcols(subject)[column] = 
              as.character(S4Vectors::mcols(subject)[[column]])
        } else if (class(S4Vectors::mcols(subject)[[column]])[1] == 
                   'DNAStringSetList') {
          S4Vectors::mcols(subject)[column] = 
            unstrsplit(CharacterList(S4Vectors::mcols(subject)[[column]]))
        }

    # Add subject metadata to query
    S4Vectors::mcols(query)[S4Vectors::queryHits(hits), 
                                paste(column, column_suffix, sep='')] = 
      S4Vectors::mcols(subject)[S4Vectors::subjectHits(hits), column]
    }
    return(query)
}

#' @export
#' @rdname overlap_variants
overlap_variants = function(object_1,
                            object_2,
                            sample_1="sample_1",
                            sample_2="sample_2") {

    # Find the intersection of all ranges in both objects
    intersect.gr = S4Vectors::intersect(object_1, object_2)

    # Add metadata from both objects to the union object
    intersect.gr = add_metadata(intersect.gr, object_1, paste0(".", sample_1))
    intersect.gr = add_metadata(intersect.gr, object_2, paste0(".", sample_2))

    # Convert to data frame
    data = GenomicRanges::as.data.frame(intersect.gr)

    # Remove non-complete variants
    alleles = paste(c("A1", "A1", "A2", "A2"), c(sample_1, sample_2), sep=".")
    data = data[complete.cases(data[, alleles]), ]

    # Add empty data frame with sample names if no variants overlap
    if (nrow(data) == 0) {
        data[1, "sample_1"] = sample_1
        data[1, "sample_2"] = sample_2
    } else {
        data$sample_1 = sample_1
        data$sample_2 = sample_2
    }
    
    # Remove redundant sample name columns
    redundant = paste0("sample.", c(sample_1, sample_2))
    data = data[!names(data) %in% redundant]

    # Return the final data frame
    return(data)
}
