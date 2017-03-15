#' Load existing packages and install missing ones.
#'
#' This is a function for automatic installation and loading of R packages
#' that are not already installed. Installation is first attempted from CRAN
#' (\url{https://cran.r-project.org/}), followed by BioConductor
#' (\url{http://bioconductor.org/}).
#'
#' @param packages The packages to be installed.
#' @examples
#' loadPackages('docopt')
#' loadPackages(c('docopt', 'plyr'))

loadPackages = function(packages) {

    # Loop over every package submitted
    for ( package in packages ) {

        # Load package, install if missing
        if ( ! require(package, character.only=TRUE) ) {

            # First try to install from CRAN ...
            tryCatch (silent=TRUE,
                install.package(package, 
                                repos='http://cran.us.r-project.org')
                
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

            # Load package after installation
            require(package)
        }
    }
}
