#' @title Read COSMIC data
#'
#' @description Read COSMIC sample-specific mutational data.
#'
#' @details This function reads the COSMIC data files (e.g.
#' "CosmicCLP_MutantExport.tsv.gz") and returns a GRanges object with all the
#' listed mutations for the specified sample, which can then be use in
#' downstream profile  comparisons. Only non-duplicated (gene-level) SNVs are
#' included in COSMIC profiles.
#'
#' @export
#' @rdname read_cosmic
#' @param file_path The COSMIC data file path (path).
#' @param sample_name Subset the data on sample name (character).
#' @param primary_site Subset the data on primary tumour site (character).
#' @return A dataframe with COSMIC SNVs.
#'
#' @examples
#' # Path to COSMIC test data
#' file <- system.file("extdata",
#'                     "subset_CosmicCLP_MutantExport.tsv.gz",
#'                     package = "seqCAT")
#'
#' # Read COSMIC test data for the HCT116 cell line
#' cosmic_hct116 <- read_cosmic(file, "HCT116")
read_cosmic <- function(file_path,
                        sample_name  = NULL,
                        primary_site = NULL) {

    # Message
    message("Reading COSMIC data ...")

    # Read COSMIC data
    cosmic <- utils::read.table(file_path,
                                header           = TRUE,
                                sep              = "\t",
                                quote            = "",
                                comment          = "",
                                fill             = TRUE,
                                stringsAsFactors = FALSE)

    # Remove sites without a listed position
    cosmic <- cosmic[cosmic$Mutation.genome.position != "", ]

    # Fix column naming
    names(cosmic) <- tolower(gsub("\\.", "_", names(cosmic)))

    # Simplify COSMIC sample/site names
    cosmic$sample_name <- toupper(gsub("[-. ]", "", cosmic$sample_name))
    cosmic$primary_site <- toupper(gsub("[-. ]", "", cosmic$primary_site))

    # Get data for selected sample (if applicable)
    if (!is.null(sample_name)) {
        sample_name <- toupper(gsub("[-. ]", "", sample_name))
        if (any(grepl(sample_name, cosmic$sample_name))) {
            cosmic <- cosmic[grep(sample_name, cosmic$sample_name), ]
            if (length(unique(cosmic$sample_name)) > 1) {
                cosmic <- cosmic[cosmic$sample_name == sample_name, ]
            }
        } else {
            stop("the sample ", sample_name, " is either not present in the data",
                 " or has no listed SNVs.")
        }
    }

    # Get data for selected tumour site (if applicable)
    if (!is.null(primary_site)) {
        primary_site <- toupper(gsub("[-. ]", "", primary_site))
        if (any(grepl(primary_site, cosmic$primary_site))) {
            cosmic <- cosmic[grep(primary_site, cosmic$primary_site), ]
            if (length(unique(cosmic$primary_site)) > 1) {
                cosmic <- cosmic[cosmic$primary_site == primary_site, ]
            }

            # Rename sample to primary_site for downstream comparisons
            cosmic$sample_name <- primary_site
        } else {
            stop("the primary site ", primary_site, " is either not present in the data",
                 " or has no listed SNVs.")
        }
    }

    # Keep only SNVs
    cosmic <- cosmic[grep("Substitution", cosmic$mutation_description), ]

    # Remove duplicate positions
    cosmic$gene_name <- gsub("_ENST\\d+", "", cosmic$gene_name)
    cosmic <- cosmic[!duplicated(cosmic[c("gene_name",
                                          "mutation_genome_position")]), ]

    # Separate chromosomes and positions
    positions <- data.frame(do.call(rbind, strsplit(as.vector(
        cosmic$mutation_genome_position), split = ":|-")))
    names(positions) <- c("chr", "start", "end")
    cosmic <- cbind(cosmic, positions)
    cosmic$start <- as.numeric(as.character(cosmic$start))
    cosmic$end <- as.numeric(as.character(cosmic$end))
    cosmic$mutation_genome_position <- NULL

    # Get REF and ALT
    mut <- as.character(cosmic$mutation_cds)
    cosmic$REF <- substr(cosmic$mutation_cds, nchar(mut) - 2, nchar(mut) - 2)
    cosmic$ALT <- substr(cosmic$mutation_cds, nchar(mut), nchar(mut))

    # Complement reverse strand mutations
    # (in order to compare with the normal all-forward VCF standard)
    idx <- row.names(cosmic[cosmic$strand == "-", ])
    cosmic[idx, "REF"] <- chartr("ATCG", "TAGC", cosmic[idx, "REF"])
    cosmic[idx, "ALT"] <- chartr("ATCG", "TAGC", cosmic[idx, "ALT"])
    cosmic$strand <- NULL

    # Get COSMIC alleles
    cosmic$A1 <- cosmic$REF
    cosmic[cosmic$mutation_zygosity == "hom", "A1"] <-
        cosmic[cosmic$mutation_zygosity == "hom", "ALT"]
    cosmic$A2 <- cosmic$ALT
    cosmic$mutation_zygosity <- NULL

    # Rename columns to adhere to non-COSMIC structure
    names(cosmic) <- gsub("start", "pos", names(cosmic))
    names(cosmic) <- gsub("gene_name", "gene", names(cosmic))
    names(cosmic) <- gsub("sample_name", "sample", names(cosmic))
    names(cosmic) <- gsub("accession_number", "ENSTID", names(cosmic))
    names(cosmic) <- gsub("mutation_", "", names(cosmic))
    cosmic$end <- NULL
    cosmic$strand <- NULL

    # Re-order columns
    first <- c("chr", "pos", "REF", "ALT", "A1", "A2")
    extra <- setdiff(names(cosmic), first)
    cosmic <- cosmic[c(first, extra)]

    # Add "COSMIC" to sample name
    cosmic$sample <- paste0("COSMIC.", cosmic$sample)

    # Return final data frame
    return(cosmic)
}

#' @title List COSMIC sample names
#'
#' @description List all available samples in the COSMIC database
#'
#' @details This function lists the available sample names in the provided
#' COSMIC file (e.g. CosmicCLP_MutantExport.tsv.gz), and takes about half the
#' time it takes to read the full file with the read_cosmic function, making it
#' useful for just seeing if your particular sample is listed in COSMIC or not.
#'
#' @export
#' @rdname list_cosmic
#' @param file_path The file containing COSMIC data (path).
#' @return A vector of sample names
#' 
#' @examples
#' file <- system.file("extdata",
#'                     "subset_CosmicCLP_MutantExport.tsv.gz",
#'                     package = "seqCAT")
#' cosmic_samples <- list_cosmic(file)
list_cosmic <- function(file_path) {

    # Message
    message("Reading COSMIC data ...")

    # Get header and the number of columns
    header <- utils::read.table(file_path, sep = "\t", nrow = 1)

    # Set colClasses to only read "Sample names" column
    col_classes <- c(rep("NULL", 4), "character",
                     rep("NULL", ncol(header) - 5))

    # Read COSMIC data
    cosmic <- utils::read.table(file_path,
                                header           = TRUE,
                                sep              = "\t",
                                quote            = "",
                                comment          = "",
                                fill             = TRUE,
                                stringsAsFactors = FALSE,
                                colClasses       = col_classes)

    # Simplify sample names
    cosmic$Sample.name <- toupper(gsub("[-. ]", "", cosmic$Sample.name))

    # Return unique sample names
    cosmic_samples <- sort(unique(cosmic$Sample.name))
    return(cosmic_samples)
}
