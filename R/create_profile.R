#' Create an SNV profile from a VCF file
#'
#' This function creates a SNV profile from a given VCF file by extracting the
#' variants that pass the filtering criterias. It can either be performed using
#' R, or by the create_profile.py function included (which requires that Python
#' is installed, along with the PyVCF package). Profile creation is performed
#' to facilitate and accelerate the cell authentication procedures, which is
#' especially relevant when more than one pairwise comparison will be performed
#' on the same sample.
#'
#' @export
#' @rdname create_profile
#' @importFrom GenomicRanges as.data.frame
#' @importFrom VariantAnnotation geno
#' @param vcf_file The VCF file from which variants will be extracted.
#' @param sample The sample in the VCF that will be extracted.
#' @param output_file Results will be output to this file
#' @param filter_depth Remove variants below this sequencing depth
#' @param python Extract variants using Python instead of R
#' @return Does not return any data object, but output results to output_file
#' @examples
#' vcf_file = system.file("extdata",
#'                        "example.vcf.gz", 
#'                        package = "CellAuthentication")
#' create_profile(vcf_file, "sample1", "profile1.txt")
#' create_profile(vcf_file, "sample1", "profile1.txt", filter_depth = 15)
#' create_profile(vcf_file, "sample1", "profile1.txt", python = TRUE)
create_profile <- function(vcf_file,
                           sample,
                           output_file,
                           filter_depth = 10,
                           python       = FALSE) {

    # Use Python
    if (python) {

        # Python script
        command <- system.file("python/create_profile.py",
                               package = "CellAuthentication")

        # Run Python code
        system2(command, args = c(vcf_file,
                                  sample,
                                  output_file,
                                  "-f", filter_depth))

    # Use R
    } else {

        # Read VCF file
        vcf <- VariantAnnotation::readVcf(vcf_file)

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
                                 remove = TRUE)

        # Add alleles
        data <- tidyr::separate_(data   = data,
                                 col    = "GT",
                                 sep    = "/",
                                 into   = c("A1", "A2"),
                                 remove = TRUE)

        data[data$A1 == 0, "A1"] <- data[data$A1 == 0, "REF"]
        data[data$A1 == 1, "A1"] <- data[data$A1 == 1, "ALT"]
        data[data$A2 == 0, "A2"] <- data[data$A2 == 0, "REF"]
        data[data$A2 == 1, "A2"] <- data[data$A2 == 1, "ALT"]

        # Initialise empty data frame for final results
        results <- data.frame(effect           = character(),
                              impact           = character(),
                              gene             = character(),
                              ENSGID           = character(),
                              feature          = character(),
                              ENSTID           = character(),
                              biotype          = character(),
                              warnings         = character(),
                              seqnames         = integer(),
                              start            = integer(),
                              rsID             = character(),
                              REF              = character(),
                              ALT              = character(),
                              DP               = integer(),
                              AD1              = integer(),
                              AD2              = integer(),
                              A1               = character(),
                              A2               = character(),
                              stringsAsFactors = FALSE)

        # Loop over each SNV
        for (n in c(1:nrow(data))) {

            # Get annotation data for current SNV
            ann <- data[n, "ANN"][[1]]

            # Separate into columns
            ann <- tidyr::separate_(as.data.frame(ann),
                                    col    = "ann",
                                    sep    = "\\|",
                                    remove = TRUE,
                                    into   = c("ALT",
                                               "effect",
                                               "impact",
                                               "gene",
                                               "ENSGID",
                                               "feature",
                                               "ENSTID",
                                               "biotype",
                                               "rank",
                                               "HGSV.c",
                                               "HGSV.p",
                                               "cDNA.pos",
                                               "CDS.pos",
                                               "protein.pos",
                                               "distance",
                                               "warnings"))

            # Remove unwanted data columns
            ann <- dplyr::select_(ann,
                                  "-ALT",
                                  "-rank",
                                  "-HGSV.c",
                                  "-HGSV.p",
                                  "-cDNA.pos",
                                  "-CDS.pos",
                                  "-protein.pos",
                                  "-distance")

            # Keep only the highest impact SNV(s)
            impacts <- unique(ann$impact)
            if ("HIGH" %in% impacts) {
                ann <- ann[ann$impact == "HIGH", ]
            } else if ("MODERATE" %in% impacts) {
                ann <- ann[ann$impact == "MODERATE", ]
            } else if ("LOW" %in% impacts) {
                ann <- ann[ann$impact == "LOW", ]
            }

            # SNV data columns
            data_cols <- c("seqnames",
                           "start",
                           "rsID",
                           "REF",
                           "ALT",
                           "DP",
                           "AD1",
                           "AD2",
                           "A1",
                           "A2")

            # Add SNV data to each annotation
            for (col in data_cols) {
                ann[[col]] <- data[n, col]
            }

            # Append to final results data frame
            results <- rbind(results, ann)

        }

        # Finalise output
        results <- results[c("seqnames",
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
                                              gsub("\\:", "\\.",
                                                   results$ENSTID),
                                              results$effect,
                                              results$feature,
                                              results$biotype), ]

        # Write results to file
        utils::write.table(results,
                           file      = output_file,
                           sep       = "\t",
                           row.names = FALSE,
                           quote     = FALSE)
    }
}
