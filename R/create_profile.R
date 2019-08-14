#' @title SNV profile creation
#'
#' @description Create an SNV profile from data in a VCF file.
#'
#' @details This function creates a SNV profile from a given VCF file by
#' extracting the variants that pass the filtering criterias. Profile creation
#' is performed to facilitate and accelerate the cell authentication
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
#' @param min_depth Filter variants below this sequencing depth (integer).
#' @param filter_vc Filter variants failing variant caller criteria (boolean).
#' @param filter_mt Filter mitochondrial variants (boolean).
#' @param filter_ns Filter non-standard chromosomes (boolean).
#' @param filter_gd Filter duplicate variants at the gene-level (boolean).
#' @param filter_pd Filter duplicate variants at the position-level (boolean).
#' @return A data frame.
#'
#' @examples
#' # Path to the test VCF file
#' vcf_file = system.file("extdata", "test.vcf.gz", package = "seqCAT")
#'
#' # Create SNV profiles
#' profile_1 <- create_profile(vcf_file, "sample1")
#' profile_1 <- create_profile(vcf_file, "sample1", min_depth = 15)
create_profile <- function(vcf_file,
                           sample,
                           min_depth = 10,
                           filter_vc = TRUE,
                           filter_mt = TRUE,
                           filter_ns = TRUE,
                           filter_gd = TRUE,
                           filter_pd = FALSE) {

    # Message
    message("Reading VCF file ...")

    # Check if provided sample is present in VCF
    vcf_header <- VariantAnnotation::scanVcfHeader(vcf_file)
    vcf_samples <- VariantAnnotation::samples(vcf_header)
    if (!(sample %in% vcf_samples)) {
        stop(paste0("Sample \"", sample, "\" is not present in the VCF file"))
    }

    # Define VCF parameters to be read
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

    # Stop execution if DP, AD or GT is empty
    for (metacol in c("DP", "AD", "GT")) {
        if (all(is.na(S4Vectors::mcols(gr)[[metacol]]))) {
            stop(paste("VCF contains no", metacol, "data"))
        }
    }

    # Add annotations, if present
    if (annotations) {
        gr$ANN <- VariantAnnotation::info(vcf)$ANN
    }

    # Set ALT as character
    gr$ALT <- S4Vectors::unstrsplit(IRanges::CharacterList(gr$ALT))

    # Filter variants
    data <- filter_variants(gr,
                            min_depth = min_depth,
                            filter_vc = filter_vc,
                            filter_mt = filter_mt,
                            filter_ns = filter_ns)

    # Get rsIDs, if present
    data$rsID <- row.names(data)
    row.names(data) <- NULL
    data[!grepl("^rs[0-9]+", data$rsID), "rsID"] <- "None"

    # Separate allelic depths
    data$AD <- gsub("c\\(", "", gsub("\\)", "", data$AD))
    data <- tidyr::separate(data   = data,
                            col    = "AD",
                            into   = c("AD1", "AD2"),
                            fill   = "right",
                            remove = TRUE)

    # Add alleles
    data <- tidyr::separate(data   = data,
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
    order <- c("chr",
               "pos",
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
               "warnings")
    order <- order[order %in% names(data)]
    data <- data[order]

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

    # Remove duplicate rows, if present
    data <- unique(data)

    data <- filter_duplicates(data,
                              filter_gd = filter_gd,
                              filter_pd = filter_pd)

    # Add sample to profile
    data$sample <- sample

    # Return profile as data frame
    data <- as.data.frame(data)
    return(data)
}

# Function for filtering VCF annotations
filter_annotations <- function(data) {

    # Separate ANN into rows
    data <- tidyr::unnest(data, !!rlang::sym("ANN"))

    # Separate ANN into columns
    data <- tidyr::separate(data,
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
    to_remove <- c("-ALT2",
                   "-rank",
                   "-HGSV_c",
                   "-HGSV_p",
                   "-cDNA_pos",
                   "-CDS_pos",
                   "-protein_pos",
                   "-distance")
    data <- data[, !(names(data) %in% to_remove)]

    # Impact factor priority
    priority <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

    # Initialise vector for storing row indexes to keep after impact filtration
    unique_pos <- unique(data$pos)
    unique_pos_len <- length(unique_pos)
    idx <- vector(mode = "list", length = unique_pos_len)

    # Loop over each position and get row indexes of lower impacts
    for (n in seq_len(unique_pos_len)) {

        # Current position
        current <- data[data$pos == unique_pos[n], ]

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

#' @title SNV profile creation
#'
#' @description Create SNV profiles from all VCF files in a directory
#'
#' @details This functions is a convenience-wrapper for the `create_profile`
#' function, which will create SNV profiles for each and every VCF file in the
#' provided directory. The file naming scheme used is `<sample>.vcf` and will
#' dictate the each profile's sample name.
#'
#' @export
#' @rdname create_profiles
#' @param vcf_dir The VCF directory from which the profiles will be created
#'  (path).
#' @param min_depth Remove variants below this sequencing depth (integer).
#' @param filter_vc Filter variants failing variant caller criteria (boolean).
#' @param filter_mt Filter mitochondrial variants (boolean).
#' @param filter_ns Filter non-standard chromosomes (boolean).
#' @param filter_gd Filter duplicate variants at the gene-level (boolean).
#' @param filter_pd Filter duplicate variants at the position-level (boolean).
#' @param pattern Only create profiles for a subset of files corresponding to
#'  this pattern (character).
#' @param recursive Find VCF files recursively in sub-directories as well
#'  (boolean).
#' @return A list of data frames.
#'
#' @examples
#' # Path to the test VCF directory
#' vcf_dir = system.file("extdata", package = "seqCAT")
#'
#' # Create SNV profiles
#' profiles <- create_profiles(vcf_dir, pattern = "test", recursive = TRUE)
create_profiles <- function(vcf_dir,
                            min_depth = 10,
                            filter_vc = TRUE,
                            filter_mt = TRUE,
                            filter_ns = TRUE,
                            filter_gd = TRUE,
                            filter_pd = FALSE,
                            pattern   = NULL,
                            recursive = FALSE) {

    # List VCF files to create profiles for
    files <- list.files(path       = vcf_dir,
                        full.names = TRUE,
                        pattern    = pattern,
                        recursive  = recursive)

    # Get only VCF files
    vcf_files <- grep(".vcf", files, value = TRUE)

    # Initialise profile list
    profile_list <- list()

    # Create SNV profiles for each VCF file
    for (vcf in vcf_files) {

        # Get current sample
        sample <- strsplit(vcf, "/")[[1]]
        len <- length(sample)
        sample <- strsplit(sample[len], "\\.")[[1]][1]

        # Create profile for current input file
        message("Creating profile for sample ", sample, " in VCF file ",
                vcf, " ...")
        tryCatch({
            profile <- suppressMessages(create_profile(vcf_file  = vcf,
                                                       sample    = sample,
                                                       min_depth = min_depth,
                                                       filter_vc = filter_vc,
                                                       filter_mt = filter_mt,
                                                       filter_ns = filter_ns,
                                                       filter_gd = filter_gd,
                                                       filter_pd = filter_pd))
        }, error = function(e) {
            message("ERROR; continuing to the next sample.")
        })

        # Append profile to profile list
        profile_list[[length(profile_list) + 1]] <- profile
    }

    # Return the final profile list
    return(profile_list)
}
