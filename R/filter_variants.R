#' @title Variant filtering
#'
#' @description Filter variants on several criteria.
#'
#' @details This is a function for filtering SNV profiles on several criteria:
#' sequencing depth, variant caller-specific filtering, mitochondrial variants
#' and variants in non-standard chromosomes. Only filters by sequencing depth
#' by default.
#'
#' @export
#' @param data The dataframe containing the variant data to be filtered.
#' @param min_depth Threshold for variant depth (integer).
#' @param filter_vc Filter variants not passing filtering criteria (boolean).
#' @param filter_mt Filter mitochondrial variants (boolean).
#' @param filter_ns Filter non-standard chromosomes (boolean).
#' @return A data frame containing the filtered variants.
#'
#' @examples
#' # Load test comparisons
#' data(test_profile_1)
#'
#' # Filter variants
#' filtered <- filter_variants(test_profile_1, min_depth = 15)
filter_variants <- function(data,
                            min_depth = 10,
                            filter_vc  = FALSE,
                            filter_mt = FALSE,
                            filter_ns = FALSE) {

    # Convert to GRanges object, if applicable
    if (is(data, "data.frame")) {
        gr <- df_to_gr(data)
    } else {
        gr <- data
    }

    # Remove variants not passing variant calling filters (if applicable)
    if (filter_vc) {

        # Stop execution if FILTER is empty
        if (length(gr) == length(gr[gr$FILTER == ".", ])) {
            stop(paste("VCF contains no FILTER data; please filter the VCF",
                       "or set `filter = FALSE` to ignore variant filtration"))
        }

        # Filter the data
        gr <- gr[gr$FILTER == "PASS", ]
    }

    # Remove variants below the given depth threshold
    gr <- gr[gr$DP >= min_depth & !is.na(gr$DP), ]
    
    # Remove "chr" from seqlevels
    GenomeInfoDb::seqlevels(gr) <- gsub("chr", "", GenomeInfoDb::seqlevels(gr))
    
    # Remove mitochondrial variants, if applicable
    if (filter_mt) {
        gr <- GenomeInfoDb::dropSeqlevels(gr, "MT", pruning.mode = "coarse")
    }

    # Remove non-standard chromosomes, if applicable
    if (filter_ns) {
        gr <- GenomeInfoDb::keepStandardChromosomes(gr,
                                                    pruning.mode = "coarse")
    }

    # Convert to data frame
    data <- GenomicRanges::as.data.frame(gr)

    # Remove unwanted columns
    to_remove <- c("end",
                   "width",
                   "strand",
                   "paramRangeID",
                   "QUAL")
    data <- data[, !(names(data) %in% to_remove)]
    names(data) <- c("chr", "pos", names(data)[3:ncol(data)])

    # Check for <NON_REF> sites (i.e. input may be a gVCF file)
    if ("<NON_REF>" %in% data$ALT) {
        
        # Check for non-<NON_REF> ALT alleles
        non_refs <- nrow(data[data$ALT == "<NON_REF>", ])
        if (nrow(data) == non_refs) {

            # Only <NON_REF> alleles: stop and issue error
            stop("VCF only contains <NON_REF> alleles; input may be a gVCF")
        
        } else {

            # Some <NON_REF> alleles: isuee warning and keep confident alleles
            warning(paste("VCF contains", non_refs, "/", nrow(data),
                          "<NON_REF> alleles; input may be a gVCF"))
            data <- data[data$ALT != "<NON_REF>", ]
            data$ALT <- gsub("<NON_REF>", "", data$ALT)
        }
    }

    # Remove non-SNVs
    data <- data[nchar(data$REF) == 1 &
                 nchar(data$ALT) == 1, ]

    # Check if profile contains a non-zero number of variants
    if (nrow(data) == 0) {
        stop(paste("No variants left after filtering with the current",
                   "criteria; please re-run with less stringent criteria or",
                   "skip the current sample"))
    }

    # Return the filtered data
    return(data)
}

#' @title Variant de-duplication
#'
#' @description Filter duplicated variants.
#'
#' @details This is a function for filtering duplicated variants either on the
#'  gene-level or the position-level.
#'
#' @export
#' @param data The dataframe containing the variant data to be filtered.
#' @param filter_gd Filter duplicate variants at the gene-level (boolean).
#' @param filter_pd Filter duplicate variants at the position-level (boolean).
#' @return A data frame containing the filtered variants.
#'
#' @examples
#' # Load test comparisons
#' data(test_profile_1)
#'
#' # Filter variants
#' filtered_gene <- filter_duplicates(test_profile_1)
#' filtered_position <- filter_duplicates(test_profile_1, filter_pd = TRUE)
filter_duplicates <- function(data,
                              filter_gd = TRUE,
                              filter_pd = FALSE) {

    # Remove duplicate variants at the gene-level, if applicable
    if (filter_gd) {
        if ("ENSGID" %in% names(data)) {
            data <- data[!duplicated(data[, c("chr", "pos", "ENSGID")]), ]
        } else {
            stop("No ENSGID gene data available for gene-level de-duplication")
        }
    }

    # Remove duplicate variants at the position-level, if applicable
    if (filter_pd) {
        data <- data[!duplicated(data[, c("chr", "pos")]), ]
    }

    # Return de-duplicated variants
    return(data)
}
