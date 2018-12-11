#' @title SNV profile creation
#'
#' @description Create an SNV profile from data in a VCF file.
#'
#' @details This function creates a SNV profile from a given VCF file by
#' extracting the variants that pass the filtering criterias. It can either be
#' performed using R, or by the create_profile.py function included (which
#' requires that Python  is installed, along with the PyVCF package). Profile
#' creation is performed to facilitate and accelerate the cell authentication
#' procedures, which is especially relevant when more than one pairwise
#' comparison will be performed on the same sample.
#'
#' @export
#' @rdname create_profile
#' @importFrom GenomicRanges as.data.frame
#' @importFrom VariantAnnotation geno
#' @param vcf_file The VCF file from which the profile will be created (path).
#' @param sample The sample in the VCF for which a profile will be created
#'  (character).
#' @param output_file The output file with the SNV profile (path).
#' @param filter_depth Remove variants below this sequencing depth (integer).
#' @param python Extract variants using Python instead of R (boolean).
#' @return Does not return any data object, but outputs results to output_file
#'  (to save computational time from having to repeatedly create profiles).
#'
#' @examples
#' # Path to the test VCF file
#' vcf_file = system.file("extdata", "test.vcf.gz", package = "seqCAT")
#'
#' # Create SNV profiles
#' \dontrun{
#'  create_profile(vcf_file, "sample1", "profile1.txt")
#'  create_profile(vcf_file, "sample1", "profile1.txt", filter_depth = 15)
#'  create_profile(vcf_file, "sample1", "profile1.txt", python = TRUE)
#' }
create_profile <- function(vcf_file,
                           sample,
                           output_file,
                           filter_depth = 10,
                           filter_vc    = TRUE,
                           python       = FALSE) {

    # Choose language
    if (python) {

        # Use Python
        create_profile_python(vcf_file, sample, output_file, filter_depth)

    } else {

        # Use R
        create_profile_R(vcf_file, sample, output_file, filter_depth, filter_vc)

    }
}

#' @title SNV profile creation
#'
#' @description Create SNV profiles from all VCF files in a directory
#'
#' @details This functions is a convenience-wrapper for the `create_profile`
#' function, which will create SNV profiles for each and every VCF file in the
#' provided directory. The file naming scheme used is `<sample>.vcf` and will
#' dictate the output profile filenames.
#'
#' @export
#' @rdname create_profiles
#' @param vcf_dir The VCF directory from which the profiles will be created 
#'  (path).
#' @param output_dir The output directory to put the SNV profiles in (path).
#' @param pattern Only create profiles for a subset of files corresponding to
#'  this pattern (character).
#' @param recursive Find VCF files recursively in sub-directories as well
#'  (boolean).
#' @param filter_depth Remove variants below this sequencing depth (integer).
#' @param python Extract variants using Python instead of R (boolean).
#' @return Does not return any data object, but outputs results to output_dir
#'  (to save computational time from having to repeatedly create profiles).
#'
#' @examples
#' # Path to the test VCF directory
#' vcf_dir = system.file("extdata", package = "seqCAT")
#'
#' # Create SNV profiles
#' \dontrun{
#'  create_profiles(vcf_dir, output_dir = "profiles")
#'  create_profiles(vcf_dir, output_dir = "profiles", pattern = "test")
#'  create_profiles(vcf_dir, output_dir = "profiles", recursive = TRUE)
#' }
create_profiles <- function(vcf_dir,
                            output_dir   = ".",
                            pattern      = NULL,
                            recursive    = FALSE,
                            filter_depth = 10,
                            filter_vc    = TRUE,
                            python       = FALSE) {

    # List VCF files to create profiles for
    files <- list.files(path       = vcf_dir,
                        full.names = TRUE,
                        pattern    = pattern,
                        recursive  = recursive)

    # Get only VCF files
    vcf_files <- grep(".vcf", files, value = TRUE)
    
    # Create SNV profiles for each VCF file
    for (vcf in vcf_files) {

        # Get current sample
        sample <- strsplit(vcf, "/")[[1]]
        len <- length(sample)
        sample <- strsplit(sample[len], "\\.")[[1]][1]

        # Set current output
        output <- paste0(output_dir, "/", sample, ".profile.txt")

        # Create profile for current input file
        message("Creating profile for sample ", sample, " in VCF file ",
                vcf, " ...")
        tryCatch({
            suppressMessages(create_profile(vcf,
                                            sample,
                                            output,
                                            filter_depth,
                                            filter_vc,
                                            python))
        }, error = function(e) {
            message("ERROR; continuing to the next sample.")
        })
    }
}

# Function for creating profiles with R
create_profile_R <- function(vcf_file,
                             sample,
                             output_file,
                             filter_depth,
                             filter_vc) {

    # Message
    message("Reading VCF file ...")

    # Define VCF parameters to be read
    vcf_header <- VariantAnnotation::scanVcfHeader(vcf_file)
    if ("ANN" %in% row.names(VariantAnnotation::info(vcf_header))) {
        annotations <- TRUE
        svp <- VariantAnnotation::ScanVcfParam(info = "ANN",
                                               geno = c("DP", "AD", "GT"))
    } else {
        annotations <- FALSE
        svp <- VariantAnnotation::ScanVcfParam(geno = c("DP", "AD", "GT"))
    }

    # Read VCF file
    vcf <- VariantAnnotation::readVcf(vcf_file, param = svp)

    # Message
    message("Creating SNV profile ...")

    # Gather relevant information to data GRanges object
    gr <- SummarizedExperiment::rowRanges(vcf)
    gr$DP <- as.data.frame(VariantAnnotation::geno(vcf)$DP)[[sample]]
    gr$AD <- as.data.frame(VariantAnnotation::geno(vcf)$AD)[[sample]]
    gr$GT <- as.data.frame(VariantAnnotation::geno(vcf)$GT)[[sample]]

    # Add annotations, if present
    if (annotations) {
        gr$ANN <- VariantAnnotation::info(vcf)$ANN
    }

    # Set ALT as character
    gr$ALT <- S4Vectors::unstrsplit(IRanges::CharacterList(gr$ALT))

    # Remove variants not passing variant calling filters (if applicable)
    if (filter_vc) {
        gr <- gr[gr$FILTER == "PASS", ]
    }
    gr$FILTER <- NULL

    # Remove variants below the given depth threshold
    gr <- gr[gr$DP >= filter_depth & !is.na(gr$DP), ]

    # Convert to data frame
    data <- GenomicRanges::as.data.frame(gr)

    # Check for <NON_REF> sites (i.e. input may be a gVCF file)
    if ("<NON_REF>" %in% data$ALT) {
        
        # Check for non-<NON_REF> ALT alleles
        non_refs <- nrow(data[data$ALT == "<NON_REF>", ])
        if (nrow(data) == non_refs) {

            # Only <NON_REF> alleles: stop and issue error
            stop("VCF only contains <NON_REF> alleles; input may be a gVCF")
        
        } else {

            # Some <NON_REF> alleles: isuee warning and keep confident alleles
            warning(paste("VCF contains", non_refs, "/", nrow(data),
                          "<NON_REF> alleles; input may be a gVCF"))
            data <- data[data$ALT != "<NON_REF>", ]
            data$ALT <- gsub("<NON_REF>", "", data$ALT)
        }
    }

    # Remove non-SNVs
    data <- data[nchar(data$REF) == 1 &
                 nchar(data$ALT) == 1, ]

    # Get rsIDs if existing
    data$rsID <- row.names(data)
    data[!grepl("^rs[0-9]+", data$rsID), "rsID"] <- "None"

    # Remove unwanted columns
    row.names(data) <- NULL
    data <- dplyr::select_(data,
                           "-end",
                           "-width",
                           "-strand",
                           "-paramRangeID",
                           "-QUAL")

    # Separate allelic depths
    data$AD <- gsub("c\\(", "", gsub("\\)", "", data$AD))
    data <- tidyr::separate_(data   = data,
                             col    = "AD",
                             into   = c("AD1", "AD2"),
                             fill   = "right",
                             remove = TRUE)

    # Add alleles
    data <- tidyr::separate_(data   = data,
                             col    = "GT",
                             sep    = "/",
                             into   = c("A1", "A2"),
                             fill   = "right",
                             remove = TRUE)

    data[data$A1 == 0, "A1"] <- data[data$A1 == 0, "REF"]
    data[data$A1 == 1, "A1"] <- data[data$A1 == 1, "ALT"]
    data[data$A2 == 0, "A2"] <- data[data$A2 == 0, "REF"]
    data[data$A2 == 1, "A2"] <- data[data$A2 == 1, "ALT"]

    # Separate and filter annotations, if present
    if (annotations) {
        data <- filter_annotations(data)
    }

    # Re-order output
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
               "warnings")
    order <- order[order %in% names(data)]
    data <- data[order]

    names(data) <- c("chr", "pos", names(data)[3:ncol(data)])

    # Remove duplicate rows, if present
    data <- unique(data)

    # Sort output
    if (annotations) {
        data <- data[order(as.character(data$chr),
                                 data$pos,
                                 data$gene,
                                 data$ENSGID,
                                 gsub("\\:", "\\.", data$ENSTID),
                                 data$effect,
                                 data$feature,
                                 data$biotype), ]
    } else {
        data <- data[order(as.character(data$chr),
                                 data$pos), ]
    }

    # Write data to file
    utils::write.table(data,
                       file      = output_file,
                       sep       = "\t",
                       row.names = FALSE,
                       quote     = FALSE)
    # Output message
    message("Created and stored SNV profile for \"", sample, "\" in [",
            output_file, "].")
}

# Function for filtering VCF annotations
filter_annotations <- function(data) {

    # Separate ANN into rows
    data <- dplyr::mutate_(data, .dots = stats::setNames(strsplit("ANN", ", "),
                                                         "ANN"))
    data <- tidyr::unnest_(data, "ANN")

    # Separate ANN into columns
    data <- tidyr::separate_(data,
                            col    = "ANN",
                            sep    = "\\|",
                            extra  = "drop",
                            fill   = "right",
                            remove = TRUE,
                            into   = c("ALT2",
                                       "effect",
                                       "impact",
                                       "gene",
                                       "ENSGID",
                                       "feature",
                                       "ENSTID",
                                       "biotype",
                                       "rank",
                                       "HGSV_c",
                                       "HGSV_p",
                                       "cDNA_pos",
                                       "CDS_pos",
                                       "protein_pos",
                                       "distance",
                                       "warnings"))

    # Remove unwanted data columns
    data <- dplyr::select_(data,
                          "-ALT2",
                          "-rank",
                          "-HGSV_c",
                          "-HGSV_p",
                          "-cDNA_pos",
                          "-CDS_pos",
                          "-protein_pos",
                          "-distance")

    # Impact factor priority
    priority <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

    # Initialise vector for storing row indexes to keep after impact filtration
    unique_pos <- unique(data$start)
    unique_pos_len <- length(unique_pos)
    idx <- vector(mode = "list", length = unique_pos_len)

    # Loop over each position and get row indexes of lower impacts
    for (n in seq_len(unique_pos_len)) {

        # Current position
        current <- data[data$start == unique_pos[n], ]

        # Skip if current position only contains a single entry
        if (nrow(current) == 1) {
            next
        }

        # Skip if current position only contains a single impact category
        if (length(unique(current$impact)) == 1) {
            next
        }

        # Find the highest impact in the current position
        highest <- min(which(priority %in% unique(current$impact) == TRUE))

        # Remove rows in data that don't contain the highest impact
        idx[[n]] <- row.names(current[current$impact != priority[highest], ])
    }

    # Remove impacts
    idx <- unlist(idx)
    data <- data[!(row.names(data) %in% idx), ]

    # Return data
    return(data)
}

# Function for creating profiles with Python
create_profile_python <- function(vcf_file,
                                  sample,
                                  output_file,
                                  filter_depth = 10) {

    # Message
    message("Creating SNV profile with Python ...")

    # Python script
    command <- system.file("python/create_profile.py", package = "seqCAT")

    # Run Python code
    out <- suppressWarnings(
        try(silent = TRUE, system2(command = command,
                                   stderr  = TRUE,
                                   args    = c(vcf_file,
                                               sample,
                                               output_file,
                                               "-f", filter_depth)))
        )

    # Catch Python errors and print to user
    if (length(out) > 0) {
        message("The Python script produced an ERROR; please make sure that ",
                "Python and the PyVCF module are installed.")
        message("Original error message:")
        write(out[seq_len(out)], file = "")
    } else {
        message("Created and stored SNV profile for \"", sample, "\" in [",
                output_file, "].")
    }
}
