#' Find overlapping variants in two GenomicRanges objects.
#'
#' This are functions for finding overlapping variants in two different
#' GenomicRanges objects, returning another GRanges object.
#'
#' The add_metadata function is a function for adding metadata (i.e. any
#' column that is not the "seqnames", "start" or "end" fields in a
#' GenomicRanges object) from the "subject" GRanges object to the "query" 
#' GRanges object. The overlap_variants function, on the other hand, is a
#' wrapper function that calls add_metadata twice in succession, ensuring that
#' even non-overlapping variants are included in the final result.
#'
#' @param query The query to add subject metadata to.
#' @param subject The subject whose metadata gets added to query.
#' @param column_suffix A string which will be added to all the metadata
#'     columns from the current add_metadata operation
#' @param object_1 The first variant GRanges object.
#' @param object_2 The second variant GRanges object.
#' @return Each function returns a single GRanges object.
#' @examples
#' add_metadata(data_first, data_second)
#' overlap_variants(data_first, data_second)

#' @export
#' @rdname overlap_variants
add_metadata = function(query, subject, column_suffix) {

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
overlap_variants = function(object_1, object_2) {

    # Find the union of all ranges in both objects
    union.gr = S4Vectors::union(object_1, object_2)

    # Add metadata from both objects to the union object
    union.gr = add_metadata(union.gr, object_1, '.input_1')
    union.gr = add_metadata(union.gr, object_2, '.input_2')

    # Convert to data frame
    data = GenomicRanges::as.data.frame(union.gr)

    # Remove non-complete genotypes
    alleles = c("A1.input_1", "A2.input_1", "A1.input_2", "A2.input_2")
    data = data[complete.cases(data[, alleles]), ]
    
    # Return the final data frame
    return(data)
}
