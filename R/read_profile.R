#' @title Read SNV profile
#'
#' @description \link{read_profile} reads SNV profiles as GRanges objects for
#'  downstream profile comparisons
#'
#' @details This is a function for reading SNV profiles extracted from VCF
#' files.  The data is returned as a GenomicRanges object, suitable for merging
#' of metadata.
#' 
#' @export
#' @param file The file path to be read.
#' @param sample_name The sample name from which the file originates.
#' @return A GenomicRanges object.
#'
#' @examples
#' # Path to test data
#' profile = system.file("extdata",
#'                       "test_profile_1.txt.gz",
#'                       package = "seqCAT")
#' 
#' # Read test profile
#' profile_1 <- read_profile(profile, "sample1")
read_profile <- function(file, sample_name) {

    # Message
    message("Reading profile for ", sample_name, " in file ", basename(file),
            " ...")

    # Read data
    data <- utils::read.table(file             = file,
                              sep              = "\t",
                              header           = TRUE,
                              stringsAsFactors = FALSE)

    # Remove duplicate variants
    data <- data[!duplicated(data[, c("chr", "pos", "ENSGID")]), ]

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
    data_gr <- GenomeInfoDb::dropSeqlevels(data_gr,
                                           "MT",
                                           pruning.mode = "coarse")
    data_gr <- GenomeInfoDb::keepStandardChromosomes(data_gr)

    # Return the GenomicRanges object
    return(data_gr)
}
