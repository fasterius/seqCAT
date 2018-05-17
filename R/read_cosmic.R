#' @title Read COSMIC SNV data
#'
#' @description Read COSMIC cell line-specific mutational data.
#'
#' @details This function reads the "CosmicCLP_MutantExport.tsv.gz" file
#' obtained from COSMIC and returns a GRanges object with all the listed
#' mutations for the specified cell line, which can then be use in downstream
#' profile  comparisons. Only non-duplicated (gene-level) SNVs are included in
#' COSMIC profiles.
#'
#' @export
#' @rdname read_cosmic
#' @param file_path The CosmicCLP_MutantExport.tsv.gz file (path).
#' @param cell_line The cell line to be investigated (character).
#' @return A GRanges object with COSMIC SNVs.
#'
#' @examples
#' # Path to COSMIC test data
#' file <- system.file("extdata",
#'                     "subset_CosmicCLP_MutantExport.tsv.gz",
#'                     package = "seqCAT")
#'
#' # Read COSMIC test data for HCT116 cell line
#' cosmic_hct116 <- read_cosmic(file, "HCT116")
read_cosmic <- function(file_path, cell_line) {

    # Message
    message("Reading COSMIC cell line data ...")

    # Read COSMIC data
    cosmic <- utils::read.table(file_path,
                                header           = TRUE,
                                sep              = "\t",
                                quote            = "\"",
                                comment          = "",
                                stringsAsFactors = FALSE)

    # Keep only relevant columns
    cosmic <- cosmic[c("Gene.name",
                       "Sample.name",
                       "Mutation.ID",
                       "Mutation.CDS",
                       "Mutation.AA",
                       "Mutation.Description",
                       "Mutation.zygosity",
                       "Mutation.genome.position",
                       "strand",
                       "Mutation.somatic.status",
                       "Mutation.verification.status")]
    names(cosmic) <- tolower(gsub("\\.", "_", names(cosmic)))

    # Simplify COSMIC cell line names
    cosmic$sample_name <- toupper(gsub("[-. ]", "", cosmic$sample_name))
    cell_line <- toupper(gsub("[-. ]", "", cell_line))

    # Check if cell line is available in the data
    if (!any(grepl(cell_line, cosmic$sample_name))) {
        stop("the cell line ", cell_line, " is not available in COSMIC.")
    }

    # Remove sites without a listed position
    cosmic <- cosmic[cosmic$mutation_genome_position != "", ]

    # Get data for selected cell line
    cosmic <- cosmic[grep(cell_line, cosmic$sample_name), ]
    if (length(unique(cosmic$sample_name)) > 1) {
        cosmic <- cosmic[cosmic$sample_name == cell_line, ]
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
    names(cosmic) <- c("gene",
                       "sample",
                       "ID",
                       "CDS",
                       "AA",
                       "description",
                       "somatic_status",
                       "verification_status",
                       "chr",
                       "start",
                       "end",
                       "REF",
                       "ALT",
                       "A1",
                       "A2")

    # Add "COSMIC" to sample name
    cosmic$sample <- paste0("COSMIC.", cosmic$sample)

    # Convert to GRanges object
    cosmic_gr <- GenomicRanges::makeGRangesFromDataFrame(
        cosmic,
        keep.extra.columns      = TRUE,
        ignore.strand           = TRUE,
        seqinfo                 = NULL,
        seqnames.field          = "chr",
        start.field             = "start",
        end.field               = "end",
        starts.in.df.are.0based = FALSE)

    # Rename chromosomes (23, 24) to (X, Y)
    GenomeInfoDb::seqlevels(cosmic_gr, pruning.mode = "coarse") <-
        c(as.character(1:22), "X", "Y")

    # Return final GRanges object
    return(cosmic_gr)
}

#' @title List COSMIC cell lines
#'
#' @description List all available cell lines in the COSMIC database
#'
#' @details This function lists the available cell lines in the provided
#' CosmicCLP_MutantExport.tsv.gz file, and takes about half the time it takes
#' to read the full file with the read_cosmic function, making it useful for
#' just seeing if your particular cell line is listed in COSMIC or not.
#'
#' @export
#' @rdname list_cosmic
#' @return A vector of cell line names
#' 
#' @examples
#' file <- system.file("extdata",
#'                     "subset_CosmicCLP_MutantExport.tsv.gz",
#'                     package = "seqCAT")
#' cell_lines <- list_cosmic(file)
list_cosmic <- function(file_path) {

    # Message
    message("Reading COSMIC cell line data ...")

    # Set colClasses to only read "Sample names" column
    col_classes <- c(rep("NULL", 4), "character", rep("NULL", 33))

    # Read COSMIC data
    cosmic <- utils::read.table(file_path,
                                header           = TRUE,
                                sep              = "\t",
                                quote            = "\"",
                                comment          = "",
                                stringsAsFactors = FALSE,
                                colClasses       = col_classes)

    # Simplify cell line names
    cosmic$Sample.name <- toupper(gsub("[-. ]", "", cosmic$Sample.name))

    # Get vector of the cell lines and return it
    cell_lines <- sort(unique(cosmic$Sample.name))
    return(cell_lines)
}
