#' Find overlapping variants in two GenomicRanges objects
#'
#' This are functions for finding overlapping variants in two different
#' GenomicRanges objects, returning another GRanges object.
#'
#' The addMetadata function is a function for adding metadata (i.e. any
#' column that is not the "seqnames", "start" or "end" fields in a
#' GenomicRanges object) from the "subject" GRanges object to the "query" 
#' GRanges object. The variantOverlaps function, on the other hand, is a
#' wrapper function that calls addMetadata twice in succession, ensuring that
#' even non-overlapping variants are included in the final result.
#' 
#' @section 
#'
#' @param query The query to add subject metadata to.
#' @param subject The subject whose metadata gets added to query.
#' @param column_suffix A string which will be added to all the metadata
#'     columns from the current addMetadata operation
#' @param object_1 The first variant GRanges object.
#' @param object_2 The second variant GRanges object.
#' @return Each function returns a single GRanges object.
#' @examples
#' addMetadata(data_first, data_second)
#' variantOverlaps(data_first, data_second)

#' @rdname addMetadata
addMetadata = function(query, subject, column_suffix) {

    # Find overlapping ranges
    hits = findOverlaps(query, subject)

    for ( column in names(mcols(subject)) ) {

        # Create empty metadata column to be filled
        mcols(query)[paste(column, column_suffix, sep='')] = NA

        # Convert DNAStringSet / DNAStringSetList columns to character vectors
        if (class(mcols(subject)[[column]])[1] == 'DNAStringSet') {
          mcols(subject)[column] = as.character(mcols(subject)[[column]])
        } else if (class(mcols(subject)[[column]])[1] == 'DNAStringSetList') {
          mcols(subject)[column] = 
            unstrsplit(CharacterList(mcols(subject)[[column]]))
        }

    # Add subject metadata to query
    mcols(query)[queryHits(hits), paste(column, column_suffix, sep='')] = 
      mcols(subject)[subjectHits(hits), column]
    }
    return(query)
}

#' @rdname variantOverlaps
variantOverlaps = function(object_1, object_2) {

    # Load packages
    source('load_packages.R', chdir=TRUE)
    loadPackages('GenomicRanges')

    # Find the union of all ranges in both objects
    union.gr = union(object_1, object_2)

    # Add metadata from both objects to the union object
    union.gr = addMetadata(union.gr, object_1, '.input_1')
    union.gr = addMetadata(union.gr, object_2, '.input_2')

    # Convert to data frame
    data = as.data.frame(union.gr)
    
    # Return the final object as a data frame
    return(data)
}
