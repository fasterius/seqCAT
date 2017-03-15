#' Installs required packages if not already installed.
#'
#' This is a function for automatic installation of any R packages that are
#' not already installed. A list of required packages are submitted to
#' \code{packages} and installed from CRAN (\url{https://cran.r-project.org/})
#' and BioConductor (\url{http://bioconductor.org/}), if not already available
#' to R.
#'
#' @param packages A list of packages (with package names as strings)
#' @examples
#' checkPackages('docopt')
#' checkPackages(c('docopt','plyr'))

checkPackages = function(packages) {

    # If any required packages are missing, install them
    if ( length(setdiff(packages, rownames(installed.packages()))) > 0 ) {

        # Output a message
        cat('installing missing packages ...\n')
        
        # First try to install from CRAN ...
        tryCatch (silent=TRUE,
            install.packages(setdiff(packages, rownames(installed.packages())),
                             repos='http://cran.us.r-project.org'),
            
            # ... then from BioConductor, if unsuccessful
            warning=function(bc) {
                source("http://bioconductor.org/biocLite.R")
                biocLite(setdiff(packages, rownames(installed.packages())))
            },
            error=function(bc) {
                source("http://bioconductor.org/biocLite.R")
                biocLite(setdiff(packages, rownames(installed.packages())))
            }
        )
    }
}

