<p align="center">
    <img src="man/figures/README_seqCAT_logo.png" width="400", alt="seqCAT"/>
</p>

## Overview

[![License: MIT][1]][2] [![Build status][3]][4] [![Coverage Status][5]][6] 

The High Throughput Sequencing Cell Authentication Toolkit (**seqCAT**) is an
R-package for authenticating, evaluating and characterisation of cells using
*single nucleotide variants* (SNVs) from sequencing data. Its input data should
be on the form of [VCF files][7], *i.e.* output from variant callers such as
the [Genome Analysis ToolKit][8] and annotated with software such as
[SnpEff][9].

## Installation

`seqCAT` is available on both [Bioconductor][10] and here on GitHub. You can
install the latest, stable version from Bioconductor like so:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("seqCAT")
```

If you are interested in the development version of `seqCAT`, you can install
it from GitHub:

```r
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

# Path to the example VCF file
vcf <- system.file("extdata", "example.vcf.gz", package = "seqCAT")

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

# Create an heatmap of comparisons and their similarity scores
plot_heatmap(comparisons[[1]])
```
<p align="center">
    <img src="man/figures/README_example_1.png", alt="Example heatmap"/>
</p>

For more detailed instructions on how to use seqCAT, please see the vignette.

## Citation

If you are using seqCAT to analyse your data, please cite the
following article:

> **seqCAT: a Bioconductor R-package for variant analysis of high throughput**
> **sequencing data**
> <br/> Fasterius E. and Al-Khalili Szigyarto C.
> <br/> *F1000Research* (2018), 7:1466
> <br/> https://f1000research.com/articles/7-1466

## License

seqCAT is released with a MIT licence. seqCAT is free software: you may
redistribute it and/or modify it under the terms of the MIT license. For more
information, please see the `LICENCE` file that comes with the seqCAT package.

[1]: https://img.shields.io/badge/License-MIT-blue.svg
[2]: https://opensource.org/licenses/MIT
[3]: https://travis-ci.org/fasterius/seqCAT.svg?branch=master
[4]: https://travis-ci.org/fasterius/seqCAT
[5]: https://coveralls.io/repos/github/fasterius/seqCAT/badge.svg?branch=master
[6]: https://coveralls.io/github/fasterius/seqCAT?branch=master

[7]: http://www.internationalgenome.org/wiki/Analysis/variant-call-format
[8]: https://software.broadinstitute.org/gatk/
[9]: http://snpeff.sourceforge.net/
[10]: https://bioconductor.org/packages/release/bioc/html/seqCAT.html
