#' @title Read SNV profile
#'
#' @description Read an SNV profile stored on disk.
#'
#' @details This is a function for reading SNV profiles created from VCF
#' files that were stored on disk.
#' 
#' @export
#' @rdname read_profile
#' @param file The SNV profile to be read (path). 
#' @param sample_name Sample name to be added; overrides profile sample if it
#'  already exists (character).
#' @return A data frame.
#'
#' @examples
#' # Path to test data
#' profile = system.file("extdata",
#'                       "test_1.profile.txt.gz",
#'                       package = "seqCAT")
#' 
#' # Read test profile
#' profile <- read_profile(profile)
read_profile <- function(file,
                         sample_name = NULL) {

    # Message
    message("Reading SNV profile in file ", basename(file), " ...")

    # Get file format
    format <- rev(strsplit(gsub("\\.gz", "", file), ".", fixed = TRUE)[[1]])[1]

    # Allowed non-txt formats
    non_text <- c("bed", "gtf", "gff", "gff2", "gff3")

    # Read SNV profile
    if (format == "txt") {

        profile <- utils::read.table(file             = file,
                                     sep              = "\t",
                                     quote            = "",
                                     header           = TRUE,
                                     stringsAsFactors = FALSE)

    } else if (format %in% non_text) {

        # Read non-txt profile
        profile_gr <- rtracklayer::import(file)

        # Convert to dataframe
        profile <- GenomicRanges::as.data.frame(profile_gr)

        # Re-order and keep relevant data
        order <- c("seqnames",
                   "start",
                   "rsID",
                   "gene",
                   "ENSGID",
                   "ENSTID",
                   "REF",
                   "ALT",
                   "impact",
                   "effect",
                   "feature",
                   "biotype",
                   "DP",
                   "AD1",
                   "AD2",
                   "A1",
                   "A2",
                   "FILTER",
                   "warnings",
                   "sample")
        order <- order[order %in% names(profile)]
        profile <- profile[order]
        if (length(names(profile)) > 2) {
            names(profile) <- c("chr", "pos", names(profile)[3:ncol(profile)])
        } else {
            names(profile) <- c("chr", "pos")
        }

    } else {

        # Stop execution for unsupported format specifications
        stop(paste0("Unsupported format specification \"", format,
                    "\"; please use txt, bed, gtf, gff, gff2 or gff3"))
    }

    # Add sample name (if applicable)
    if (!is.null(sample_name)) {
        profile$sample <- sample_name
    }

    # Return the SNV profile data frame
    return(profile)
}

#' @title Read SNV profiles
#'
#' @description Read SNV profiles in a directory. 
#'
#' @details This is a wrapper function for reading multiple SNV profiles
#'  present in a directory (and its sub-directories in recursive mode).
#' 
#' @export
#' @rdname read_profiles
#' @param profile_dir The directory containing the profiles to be read (path). 
#' @param pattern Pattern for file name or extension to be read (character).
#' @param sample_names Add sample name based on file name; overrides profile
#'  sample if it already exists (boolean).
#' @return A list of data frames.
#'
#' @examples
#' # Path to test data
#' profile_dir = system.file("extdata", package = "seqCAT")
#' 
#' # Read test profiles
#' profile_list <- read_profiles(profile_dir, pattern = "profile.txt")
read_profiles <- function(profile_dir,
                          pattern      = ".profile.txt",
                          sample_names = FALSE) {

    # Get all SNV profiles in the input directory
    profiles <- list.files(profile_dir,
                           full.names = TRUE,
                           pattern    = pattern)

    # Calculate total number of profiles to read and initialise counter
    nn_tot <- length(profiles)
    nn <- 1

    # Read all profiles
    profile_list <- list()
    for (profile in profiles) {

        # Add or override sample name (if applicable)
        if (sample_names) {
            current_sample <- gsub(pattern, "", basename(profile))
        } else {
            current_sample <- NULL
        }

        # Read current profile
        message(paste0("Reading profile in file ", basename(profile),
                       " [", nn, " / ", nn_tot, "]"))
        current_profile <-
            suppressMessages(read_profile(profile,
                                          sample_name = current_sample))

        # Append profile to profile list
        profile_list[[length(profile_list) + 1]] <- current_profile

        # Increment counter
        nn <- nn + 1
    }

    # Return profile list
    return(profile_list)
}
