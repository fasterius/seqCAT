#' @title Write SNV profile
#'
#' @description Write an SNV profile to a file for later re-use. 
#'
#' @details This is a function for writing SNV profiles (created from VCF
#' files) to disk for later re-use. Several formats are allowed, including BED,
#' GTF, GFF and normal text files, which are automatically recognised based on
#' the supplied filename.
#' 
#' @export
#' @rdname write_profile
#' @param profile The SNV profile to be written (data frame). 
#' @param file The file to write to (path). 
#' @return None; writes to disk only.
#'
#' @examples
#' # Load test profile
#' data(test_profile_1)
#'
#' # Write test profile to file
#' write_profile(test_profile_1, "test_profile_1.txt")
write_profile <- function(profile,
                          file) {

    # Get format specification from file name
    format <- rev(strsplit(file, ".", fixed = TRUE)[[1]])[1]

    # Allowed non-txt formats
    non_text <- c("bed", "gtf", "gff", "gff2", "gff3")

    # Check format and write to file
    if (format == "txt") {

        # Write to text
        utils::write.table(profile,
                           file      = file,
                           sep       = "\t",
                           row.names = FALSE,
                           quote     = FALSE)

    } else if (format %in% non_text) {

        # Convert to GRanges
        profile_gr <- convert_to_gr(profile)

        # Write to specified format
        rtracklayer::export(profile_gr, file)

    } else {

        # Stop execution for unsupported format specifications
        stop(paste0("Unsupported format specification \"", format,
                    "\"; please use txt, bed, gtf, gff, gff2 or gff3"))
    }

    # Output message
    message("Stored SNV profile in ", file, ".")
}
