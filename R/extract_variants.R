#' Extract relevant variant information from a VCF file
#'
#' This function is a wrapper to the "extract_variants.py" Python script, which
#' extracts relevant variant information from a VCF file. Extraction is
#' performed to facilitate and accelerate the cell authentication procedures,
#' which is especially relevant when more than one pairwise comparison will be
#' performed on the same sample.
#'
#' @param vcf_file The VCF file from which variants will be extracted.
#' @param sample The sample in the VCF that will be extracted.
#' @param output_file Results will be output to this file
#' @return Does not return any data object, but output results to output_file
#' @examples
#' extract_variants_file(vcf_file, sample, output_file, filter_depth=10,
#'                       python=FALSE)

#' @export
#' @rdname extract_variants
extract_variants = function(vcf_file,
                            sample,
                            output_file,
                            filter_depth=10,
                            python=FALSE) {
    
    # Extract variants
    if (python) {

        # Use Python
        message('[extracting with Python]')

        # Source the Python script
        command = system.file('python/extract_variants.py',
							  package='CellAuthentication')

        # Run Python code
        system2(command, args=c(vcf_file, sample, output_file, '-f',
				filter_depth))

    } else {

        # Use R
        message('[extracting with R]')

        # Read VCF file
        vcf = VariantAnnotation::readVcf(vcf_file)

        # Gather relevant information to data GRanges object
        gr = SummarizedExperiment::rowRanges(vcf)
        gr$ANN = VariantAnnotation::info(vcf)$ANN
        gr$DP = as.data.frame(VariantAnnotation::geno(vcf)$DP)[[sample]]
        gr$AD = as.data.frame(VariantAnnotation::geno(vcf)$AD)[[sample]]
        gr$GT = as.data.frame(VariantAnnotation::geno(vcf)$GT)[[sample]]

        # Set ALT as character
        gr$ALT = S4Vectors::unstrsplit(IRanges::CharacterList(gr$ALT))

        # Remove variants not passing variant calling filters
        gr = gr[gr$FILTER == 'PASS', ]
        gr$FILTER = NULL

        # Remove variants below the given depth threshold
        gr = gr[gr$DP >= filter_depth & !is.na(gr$DP), ]

        # Convert to data frame
        data = GenomicRanges::as.data.frame(gr)

        # Remove non-SNVs
        data = data[nchar(data$REF) == 1 &
                    nchar(data$ALT) == 1, ]

        # Get rsIDs if existing
        data$rsID = row.names(data)
        data[!grepl('^rs[0-9]+', data$rsID), 'rsID'] = 'None'

        # Remove unwanted columns
        row.names(data) = NULL
        data$end = NULL
        data$width = NULL
        data$strand = NULL
        data$paramRangeID = NULL
        data$QUAL = NULL

        # Separate allelic depths
        data = tidyr::separate(data=data, col=AD, sep='\\,\\ ',
                               into=c('AD1', 'AD2'), remove=TRUE)
		data$AD1 = gsub('c\\(', '', data$AD1)
		data$AD2 = gsub('\\)', '', data$AD2)

        # Add alleles
        data = tidyr::separate(data=data, col=GT, sep='/', into=c('A1', 'A2'),
                        remove=TRUE)
        data[data$A1 == 0, 'A1'] = data[data$A1 == 0, 'REF']
        data[data$A1 == 1, 'A1'] = data[data$A1 == 1, 'ALT']
        data[data$A2 == 0, 'A2'] = data[data$A2 == 0, 'REF']
        data[data$A2 == 1, 'A2'] = data[data$A2 == 1, 'ALT']

		# Initialise empty data frame for final results
		results = data.frame(effect=character(),
			impact=character(),
			gene=character(),
			ENSGID=character(),
			feature=character(),
			ENSTID=character(),
			biotype=character(),
			warnings=character(),
			seqnames=integer(),
			start=integer(),
			rsID=character(),
			REF=character(),
			ALT=character(),
			DP=integer(),
			AD1=integer(),
			AD2=integer(),
			A1=character(),
			A2=character(), stringsAsFactors=FALSE)

		# Loop over each SNV
		for (n in c(1:nrow(data))) {

			# Get annotation data for current SNV
			ann = data[n, 'ANN'][[1]]

			# Separate into columns
			ann = tidyr::separate(as.data.frame(ann), col=ann, sep='\\|',
				into=c('ALT', 'effect', 'impact', 'gene', 'ENSGID', 'feature',
					   'ENSTID', 'biotype', 'rank', 'HGSV.c', 'HGSV.p',
					   'cDNA.pos', 'CDS.pos', 'protein.pos', 'distance',
					   'warnings'), remove=TRUE)

			# Remove unwanted data columns
			ann$ALT = NULL
			ann$rank = NULL
			ann$HGSV.c = NULL
			ann$HGSV.p = NULL
			ann$cDNA.pos = NULL
			ann$CDS.pos = NULL
			ann$protein.pos = NULL
			ann$distance = NULL

			# Keep only the highest impact SNV(s)
			impacts = unique(ann$impact)
			if ('HIGH' %in% impacts) {

				ann = ann[ann$impact == 'HIGH', ]

			} else if ('MODERATE' %in% impacts) {

				ann = ann[ann$impact == 'MODERATE', ]

			} else if ('LOW' %in% impacts) {

				ann = ann[ann$impact == 'LOW', ]

			}

			# SNV data columns
			data.cols = c('seqnames', 'start', 'rsID', 'REF', 'ALT', 'DP',
						  'AD1', 'AD2', 'A1', 'A2')

			# Add SNV data to each annotation
			for (col in data.cols) {
				ann[[col]] = data[n, col]
			}

			# Append to final results data frame
			results = rbind(results, ann)

		}

		# Finalise output
		results = results[c('seqnames', 'start', 'rsID', 'gene', 'ENSGID', 
                            'ENSTID', 'REF', 'ALT', 'impact', 'effect',
							'feature', 'biotype', 'DP', 'AD1', 'AD2', 'A1',
							'A2', 'warnings')]

		names(results) = c('chr', 'pos', names(results)[3:18])

        # Remove duplicate rows (if present)
        results = unique(results)

        # Sort output
        results = results[order(results$chr, results$pos,
                                results$ENSGID, results$ENSTID), ]

		# Write results to file
		write.table(results, output_file, sep='\t', row.names=FALSE, 
                    quote=FALSE)
    }
}
