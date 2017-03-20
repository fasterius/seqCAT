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
#' @param output_file File to store results in, when not outputting to STDOUT.
#' @return Returns the extracted variants to STDOUT or to a file. 
#' @examples
#' extract_variants(vcf_file, sample)
#' extract_variants_file(vcf_file, sample, output_file)

#' @export
#' @rdname extract_variants
extract_variants = function(vcf_file, sample) {
   sdfadfk{{{ 
    # Source the Python script
    source(system.file('python/extract_variants.py', package='CellAuthentication'))
}

command = '/Users/erikfasterius/local/scripts/authentication/inst/python/extract_variants.py'
args = c('/Users/erikfasterius/local/data/scrna/data/GSE75688/SRR2973353/03_variant_calling/SRR2973353.vcf', 'SRR2973353', 'test.output.txt')
system2(command, args=args)
