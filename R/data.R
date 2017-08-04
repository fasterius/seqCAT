#' Example variants from "sample1"
#'
#' Variants from "sample1" read with the read_variants function, for use in
#' examples and tests.
#'
#' @docType data
#' @usage data(variants_1)
#' @format A data frame with 383 rows and 20 columns:
#' \describe{
#'     \item{seqnames}{chromosome}
#'     \item{start}{start position}
#'     \item{end}{end position}
#'     \item{rsID}{mutation ID, if available}
#'     \item{gene}{associated gene}
#'     \item{ENSGID}{ensembl gene ID}
#'     \item{ENSTID}{ensembl transcript ID}
#'     \item{REF}{reference allele}
#'     \item{ALT}{alternative allele}
#'     \item{impact}{putative variant impact}
#'     \item{effect}{variant effect}
#'     \item{feature}{transcript feature}
#'     \item{biotype}{transcript biotype}
#'     \item{DP}{total variant depth}
#'     \item{AD1}{allelic depth, allele 1}
#'     \item{AD2}{allelic depth, allele 2}
#'     \item{A1}{allele 1}
#'     \item{A2}{allele 2}
#'     \item{warnings}{warnings from variant calling}
#'     \item{sample}{sample name}
#' }
"variants_1"

#' Example variants from "sample2"
#'
#' Variants from "sample2" read with the read_variants function, for use in
#' examples and tests.
#'
#' @docType data
#' @usage data(variants_2)
#' @format A data frame with 382 rows and 20 columns:
#' \describe{
#'     \item{seqnames}{chromosome}
#'     \item{start}{start position}
#'     \item{end}{end position}
#'     \item{rsID}{mutation ID, if available}
#'     \item{gene}{associated gene}
#'     \item{ENSGID}{ensembl gene ID}
#'     \item{ENSTID}{ensembl transcript ID}
#'     \item{REF}{reference allele}
#'     \item{ALT}{alternative allele}
#'     \item{impact}{putative variant impact}
#'     \item{effect}{variant effect}
#'     \item{feature}{transcript feature}
#'     \item{biotype}{transcript biotype}
#'     \item{DP}{total variant depth}
#'     \item{AD1}{allelic depth, allele 1}
#'     \item{AD2}{allelic depth, allele 2}
#'     \item{A1}{allele 1}
#'     \item{A2}{allele 2}
#'     \item{warnings}{warnings from variant calling}
#'     \item{sample}{sample name}
#' }
"variants_2"

#' Overlapping example variants
#'
#' Overlapping variants from "sample1" and "sample2" originating from the
#' example.vcf file included in the inst/extdata directory, for use in examples
#' and tests.
#'
#' @docType data
#' @usage data(overlaps)
#' @format A data frame with 51 rows and 39 columns:
#' \describe{
#'     \item{seqnames}{chromosome}
#'     \item{start}{start position}
#'     \item{end}{end position}
#'     \item{width}{variant width}
#'     \item{strand}{variant strand}
#'     \item{rsID.sample_1}{mutation ID, sample 1}
#'     \item{gene.sample_1}{associated gene, sample 1}
#'     \item{ENSGID.sample_1}{ensembl gene ID, sample 1}
#'     \item{ENSTID.sample_1}{ensembl transcript ID, sample 1}
#'     \item{REF.sample_1}{reference allele, sample 1}
#'     \item{ALT.sample_1}{alternative allele, sample 1}
#'     \item{impact.sample_1}{putative variant impact, sample 1}
#'     \item{effect.sample_1}{variant effect, sample 1}
#'     \item{feature.sample_1}{transcript feature, sample 1}
#'     \item{biotype.sample_1}{transcript biotype, sample 1}
#'     \item{DP.sample_1}{total variant depth, sample 1}
#'     \item{AD1.sample_1}{allelic depth, allele 1, sample 1}
#'     \item{AD2.sample_1}{allelic depth, allele 2, sample 1}
#'     \item{A1.sample_1}{allele 1, sample 1}
#'     \item{A2.sample_1}{allele 2, sample 1}
#'     \item{warnings.sample_1}{warnings from variant calling, sample 1}
#'     \item{rsID.sample_2}{mutation ID, sample 2}
#'     \item{gene.sample_2}{associated gene, sample 2}
#'     \item{ENSGID.sample_2}{ensembl gene ID, sample 2}
#'     \item{ENSTID.sample_2}{ensembl transcript ID, sample 2}
#'     \item{REF.sample_2}{reference allele, sample 2}
#'     \item{ALT.sample_2}{alternative allele, sample 2}
#'     \item{impact.sample_2}{putative variant impact, sample 2}
#'     \item{effect.sample_2}{variant effect, sample 2}
#'     \item{feature.sample_2}{transcript feature, sample 2}
#'     \item{biotype.sample_2}{transcript biotype, sample 2}
#'     \item{DP.sample_2}{total variant depth, sample 2}
#'     \item{AD1.sample_2}{allelic depth, allele 1, sample 2}
#'     \item{AD2.sample_2}{allelic depth, allele 2, sample 2}
#'     \item{A1.sample_2}{allele 1, sample 2}
#'     \item{A2.sample_2}{allele 2, sample 2}
#'     \item{warnings.sample_2}{warnings from variant calling, sample 2}
#' }
"overlaps"
