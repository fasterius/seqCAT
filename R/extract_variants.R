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
#' extract_variants_file(vcf_file, sample, output_file)

#' @export
#' @rdname extract_variants
extract_variants = function(vcf_file, sample, output_file) {
    
    # Source the Python script
    command = system.file('python/extract_variants.py', package='CellAuthentication')

    # Run Python code
    system2(command, args=c(vcf_file, sample, output_file))

}
