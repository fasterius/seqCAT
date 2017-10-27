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
                           python       = FALSE) {

    # Choose language
    if (python) {

        # Use Python
        create_profile_python(vcf_file, sample, output_file, filter_depth)

    } else {

        # Use R
        create_profile_R(vcf_file, sample, output_file, filter_depth)

    }
}

# Function for creating profiles with R
create_profile_R <- function(vcf_file,
                             sample,
                             output_file,
                             filter_depth) {

    # Message
    message("Reading VCF file ...")

    # Define VCF parameters to be read
    svp <- VariantAnnotation::ScanVcfParam(info = "ANN",
                                           geno = c("DP", "AD", "GT"))

    # Read VCF file
    vcf <- VariantAnnotation::readVcf(vcf_file, param = svp)

    # Message
    message("Creating SNV profile ...")

    # Gather relevant information to data GRanges object
    gr <- SummarizedExperiment::rowRanges(vcf)
    gr$ANN <- VariantAnnotation::info(vcf)$ANN
    gr$DP <- as.data.frame(VariantAnnotation::geno(vcf)$DP)[[sample]]
    gr$AD <- as.data.frame(VariantAnnotation::geno(vcf)$AD)[[sample]]
    gr$GT <- as.data.frame(VariantAnnotation::geno(vcf)$GT)[[sample]]

    # Set ALT as character
    gr$ALT <- S4Vectors::unstrsplit(IRanges::CharacterList(gr$ALT))

    # Remove variants not passing variant calling filters
    gr <- gr[gr$FILTER == "PASS", ]
    gr$FILTER <- NULL

    # Remove variants below the given depth threshold
    gr <- gr[gr$DP >= filter_depth & !is.na(gr$DP), ]

    # Convert to data frame
    data <- GenomicRanges::as.data.frame(gr)

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

    # Re-order output
    results <- data[c("seqnames",
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
                      "warnings")]

    names(results) <- c("chr", "pos", names(results)[3:18])

    # Remove duplicate rows (if present)
    results <- unique(results)

    # Sort output
    results <- results[order(as.character(results$chr),
                                          results$pos,
                                          results$gene,
                                          results$ENSGID,
                                          gsub("\\:", "\\.", results$ENSTID),
                                          results$effect,
                                          results$feature,
                                          results$biotype), ]

    # Write results to file
    utils::write.table(results,
                       file      = output_file,
                       sep       = "\t",
                       row.names = FALSE,
                       quote     = FALSE)
    # Output message
    message("Created and stored SNV profile for \"", sample, "\" in [",
            output_file, "].")
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
        write(out[1:length(out)], file = "")
    } else {
        message("Created and stored SNV profile for \"", sample, "\" in [",
                output_file, "].")
    }
}
