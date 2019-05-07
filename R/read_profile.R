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
read_profile <- function(file) {

    # Message
    message("Reading SNV profile in file ", basename(file), " ...")

    # Read SNV profile
    profile <- utils::read.table(file             = file,
                                 sep              = "\t",
                                 quote            = "",
                                 header           = TRUE,
                                 stringsAsFactors = FALSE)

    # Return the SNV profile data frame
    return(profile)
}

#' @title Read SNV profiles
#'
#' @description Read SNV profiles in a directory. 
#'
#' @details This is a wrapper function for reading multiple SNV profiles
#'  present in a directory (and its sub-directories in recursive mode). The
#'  format used is `<sample>.profile.txt`.
#' 
#' @export
#' @rdname read_profiles
#' @param profile_dir The directory containing the profiles to be read (path). 
#' @return A list of data frames.
#'
#' @examples
#' # Path to test data
#' profile_dir = system.file("extdata",
#'                           package = "seqCAT")
#' 
#' # Read test profiles
#' profile_list <- read_profiles(profile_dir)
read_profiles <- function(profile_dir) {

    # Get all SNV profiles in the input directory
    profiles <- list.files(profile_dir, full.names = TRUE,
                           pattern = "profile.txt")

    # Calculate total number of profiles to read and initialise counter
    nn_tot <- length(profiles)
    nn <- 1

    # Read all profiles
    profile_list <- list()
    for (profile in profiles) {

        # Read current profile
        message(paste0("Reading profile in file ", basename(profile),
                       " [", nn, " / ", nn_tot, "]"))
        current_profile <- suppressMessages(read_profile(profile))

        # Append profile to profile list
        profile_list[[length(profile_list) + 1]] <- current_profile

        # Increment counter
        nn <- nn + 1
    }

    # Return profile list
    return(profile_list)
}
