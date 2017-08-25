#' SNV profile 1
#'
#' SNV profile in GRanges format from "sample1", originating from the
#' test_profile_1.txt in the inst/extdata directory, for use in unit tests.
#'
#' @docType data
#' @usage data(test_profile_1)
#' @format A GRanges object with 383 elements and 17 metadata columns:
#' \describe{
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
"test_profile_1"

#' SNV profile 2
#'
#' SNV profile in GRanges format from "sample2", originating from the
#' test_profile_2.txt in the inst/extdata directory, for use in unit tests.
#'
#' @docType data
#' @usage data(test_profile_2)
#' @format A GRanges object with 382 elements and 17 metadata columns:
#' \describe{
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
"test_profile_2"

#' Overlapping and compared SNVs
#'
#' Overlapping and compared variants from "sample1" and "sample2" originating
#' from the example.vcf file included in the inst/extdata directory, for use in
#' unit tests.
#'
#' @docType data
#' @usage data(test_comparison)
#' @format A dataframe with 51 rows and 39 columns:
#' \describe{
#'     \item{chr}{chromosome}
#'     \item{pos}{SNV position}
#'     \item{DP.sample_1}{total variant depth, sample 1}
#'     \item{AD1.sample_1}{allelic depth, allele 1, sample 1}
#'     \item{AD2.sample_1}{allelic depth, allele 2, sample 1}
#'     \item{A1.sample_1}{allele 1, sample 1}
#'     \item{A2.sample_1}{allele 2, sample 1}
#'     \item{warnings.sample_1}{warnings from variant calling, sample 1}
#'     \item{DP.sample_2}{total variant depth, sample 2}
#'     \item{AD1.sample_2}{allelic depth, allele 1, sample 2}
#'     \item{AD2.sample_2}{allelic depth, allele 2, sample 2}
#'     \item{A1.sample_2}{allele 1, sample 2}
#'     \item{A2.sample_2}{allele 2, sample 2}
#'     \item{warnings.sample_2}{warnings from variant calling, sample 2}
#'     \item{sample_1}{name, sample 1}
#'     \item{sample_2}{name, sample 2}
#'     \item{match}{status of genotype comparison}
#'     \item{rsID}{mutation ID}
#'     \item{gene}{associated gene}
#'     \item{ENSGID}{ensembl gene ID}
#'     \item{ENSTID}{ensembl transcript ID}
#'     \item{REF}{reference allele}
#'     \item{ALT}{alternative allele}
#'     \item{impact}{putative variant impact}
#'     \item{effect}{variant effect}
#'     \item{feature}{transcript feature}
#'     \item{biotype}{transcript biotype}
#' }
"test_comparison"

#' Collated similarities object
#'
#' Collated similarities of multiple sample comparisons from "sample1" and
#' "sample" from the example.vcf file, for use in unit tests.
#'
#' @docType data
#' @usage data(test_similarities)
#' @format A dataframe with 3 rows and 6 columns:
#' \describe{
#'      \item{sample_1}{name of sample 1}
#'      \item{sample_2}{name of sample 2}
#'      \item{overlaps}{the number of overlaps for the comparison}
#'      \item{matches}{the number of matches for the comparison}
#'      \item{concordance}{the concordance of the profiles}
#'      \item{similarity_score}{the similarity score of the profiles}
#' }
"test_similarities"
