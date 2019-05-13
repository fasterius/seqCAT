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

#' @title Write SNV profiles
#'
#' @description Write several SNV profiles to file for later re-use.
#'
#' @details This is a wrapper function for writing multiple SNV profiles
#'  present in a directory (and its sub-directories in recursive mode).
#'
#' @export
#' @rdname write_profiles
#' @param profile_list The SNV profiles to be written (list).
#' @param format The desired file format (character).
#' @param directory The directory to write to (path).
#' @return None; writes to disk only.
#'
#' @examples
#' # Load test profiles
#' data(test_profile_1)
#' data(test_profile_2)
#' profiles <- list(test_profile_1, test_profile_2)
#'
#' # Write test profile to file
#' write_profiles(profiles, format = "TXT", directory = "./")
write_profiles <- function(profile_list,
                           format    = "TXT",
                           directory = "./") {

    # Get format
    format <- paste0(".", tolower(format))

    # Calculate total number of profiles to read and initialise counter
    nn_tot <- length(profile_list)
    nn <- 1

    # Loop through profiles
    for (profile in profile_list) {

        # Get sample for current profile
        sample <- unique(profile$sample)

        # Create output file path
        file <- paste0(directory, sample, format)

        # Write current profile
        message(paste0("Writing ", file, " [", nn, " / ", nn_tot, "]"))
        suppressMessages(write_profile(profile, file))

        # Increment counter
        nn <- nn + 1
    }
}
