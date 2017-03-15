#' Installs missing package and loads them.
#'
#' This is a function for automatic installation and loading of any R package
#' that are not already installed. Installation is first attempted from CRAN
#' (\url{https://cran.r-project.org/}), followed by BioConductor
#' (\url{http://bioconductor.org/}).
#'
#' @param package The package to be installed.
#' @examples
#' installMissingPackage('docopt')

installMissingPackage = function(package) {

    # First try to install from CRAN ...
    tryCatch (silent=TRUE,
        install.package(package, repos='http://cran.us.r-project.org'),
        
        # ... then from BioConductor, if unsuccessful
        warning = function(bc) {
            source("http://bioconductor.org/biocLite.R")
            biocLite(package)
        },
        error = function(bc) {
            source("http://bioconductor.org/biocLite.R")
            biocLite(package)
        }
    )
    
    # Load package
    require(package)
}

