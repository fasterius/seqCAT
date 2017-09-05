# seqCAT

## Overview

The High Throughput Sequencing Cell Authentication Toolkit (**seqCAT**) is an
R-package for authenticating, evaluating and characterisation of cells using
*single nucleotide variants* (SNVs)sequencing data. Its input data should be on
the form of [VCF files][1], *i.e.* output from variant callers such as the
[Genome Analysis ToolKit][2] and annotated with software such as [SnpEff][3].

[1]: http://www.internationalgenome.org/wiki/Analysis/variant-call-format
[2]: https://software.broadinstitute.org/gatk/
[3]: http://snpeff.sourceforge.net/

## Installation

```r
# Install the latest development version from GitHub
# install.packages("devtools")
devtools::install_github("fasterius/seqCAT")
```

## Usage

The general workflow of seqCAT consists of three steps:

    1.  Creation of SNV profiles
    2.  Comparisons of SNV profiles
    3.  Evaluation of profile comparisons

```r
# Load the package
library("seqCAT")

# List the example VCF file
vcf <- system.file("extdata", "example.vcf.gz",
                   package = "seqCAT")

# Create SNV profiles
create_profile(vcf, "HCT116", "hct116_profile.txt")
create_profile(vcf, "HKE3", "hke3_profile.txt")
create_profile(vcf, "RKO", "rko_profile.txt")

# Read the created profiles
hct116 <- read_profile("hct116_profile.txt", "HCT116")
hke3 <- read_profile("hke3_profile.txt", "HKE3")
rko <- read_profile("rko_profile.txt", "RKO")

# Compare all profiles to each other
comparisons <- compare_many(list(hct116, hke3, rko))

# Create an heatmap of comparisons
plot_heatmap(comparisons[[1]])
```
<p align="center">
 <img src="man/figures/README_example_1.png", alt="Example heatmap"/>
</p>

For more detailed demonstrations on how to use seqCAT, please see the vignette.

## Citation

If you are using seqCAT to analyse your samples, please cite the following
publication: 

> Fasterius, E., Raso, C., Kennedy, S., Kolch, W., Al-Khalili C. et al. (2017)
> <br/> *A novel RNA sequencing data analysis method for cell line
> authentication* <br/> PloS One, 12(2), e0171435. <br/>
> http://doi.org/10.1371/journal.pone.0171435

## License

seqCAT is released with a MIT licence. seqCAT is free software: you may
redistribute it and/or modify it under the terms of the MIT license. For more
information, please see the `LICENCE` file that comes with the seqCAT package.
