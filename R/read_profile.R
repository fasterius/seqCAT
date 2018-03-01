#' @title Read SNV profile
#'
#' @description Read SNV profiles for use in downstream comparisons.
#'
#' @details This is a function for reading SNV profiles created from VCF
#' files.  The data is returned as a GenomicRanges object, suitable for merging
#' of metadata.
#' 
#' @export
#' @param file The SNV profile to be read (path). 
#' @param sample_name The sample of the SNV profile (character).
#' @param remove_mt Remove or keep mitochondrial variants (boolean).
#' @return A GRanges object.
#'
#' @examples
#' # Path to test data
#' profile = system.file("extdata",
#'                       "test_profile_1.txt.gz",
#'                       package = "seqCAT")
#' 
#' # Read test profile
#' profile_1 <- read_profile(profile, "sample1")
read_profile <- function(file, sample_name, remove_mt = TRUE) {

    # Message
    message("Reading profile for ", sample_name, " in file ", basename(file),
            " ...")

    # Read data
    data <- utils::read.table(file             = file,
                              sep              = "\t",
                              quote            = "",
                              header           = TRUE,
                              stringsAsFactors = FALSE)

    # Remove duplicate variants, if present
    if ("ENSGID" %in% names(data)) {
        data <- data[!duplicated(data[, c("chr", "pos", "ENSGID")]), ]
    } else {
        data <- data[!duplicated(data[, c("chr", "pos")]), ]
    }

    # Add sample name
    data$sample <- sample_name

    # Convert to GRanges object
    data_gr <- GenomicRanges::makeGRangesFromDataFrame(data,
        keep.extra.columns      = TRUE,
        ignore.strand           = TRUE,
        seqinfo                 = NULL,
        seqnames.field          = "chr",
        start.field             = "pos",
        end.field               = "pos",
        starts.in.df.are.0based = FALSE)

    # Rename and remove seqlevels
    GenomeInfoDb::seqlevels(data_gr) <-
        gsub("chr", "", GenomeInfoDb::seqlevels(data_gr))
    data_gr <- GenomeInfoDb::keepStandardChromosomes(data_gr,
                                                     pruning.mode = "coarse")

    # Remove mitochondrial variants, if applicable
    if (remove_mt) {
        data_gr <- GenomeInfoDb::dropSeqlevels(data_gr,
                                               "MT",
                                               pruning.mode = "coarse")
    }

    # Check if profile contains no variants after filtering
    if (length(data_gr) == 0) {

        # Add dummy variant to contain sample name
        dummy_variant <- GenomicRanges::GRanges("1",
                                                IRanges::IRanges(start = 0,
                                                                 end   = 0))
        S4Vectors::mcols(dummy_variant)["sample"] <- sample_name
        S4Vectors::mcols(dummy_variant)["A1"] <- NA
        S4Vectors::mcols(dummy_variant)["A2"] <- NA

        # Add dummy variant to empty profile
        data_gr <- append(data_gr, dummy_variant)
    }

    # Return the GenomicRanges object
    return(data_gr)
}
