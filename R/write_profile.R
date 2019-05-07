#' @title Write SNV profile
#'
#' @description Write an SNV profile to a file for later re-use. 
#'
#' @details This is a function for writing SNV profiles (created from VCF
#' files) to disk for later re-use.
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
    
    # Write to file
    utils::write.table(profile,
                       file      = file,
                       sep       = "\t",
                       row.names = FALSE,
                       quote     = FALSE)

    # Output message
    message("Stored SNV profile in ", file, ".")
}
