<p align="center">
    <img src="man/figures/README_seqCAT_logo.png" width="400", alt="seqCAT"/>
</p>

## Overview

[![Anaconda Cloud version][1]][2] [![License: MIT][3]][4] [![Build status][5]][6] [![Coverage Status][7]][8]

The High Throughput Sequencing Cell Authentication Toolkit (**seqCAT**) is an
R-package for authenticating, evaluating and characterisation of cells using
*single nucleotide variants* (SNVs) from sequencing data. Its input data should
be on the form of [VCF files][9], *i.e.* output from variant callers such as
the [Genome Analysis ToolKit][10] and annotated with software such as
[SnpEff][11].

## Installation

The `seqCAT` package is available on both [Bioconductor][12] and here on
GitHub. You can install the latest, stable version from Bioconductor like so:

```r
# install.packages("BiocManager")
BiocManager::install("seqCAT")
```

If you are interested in the development version of `seqCAT`, you can install
it from GitHub:

```r
# install.packages("devtools")
devtools::install_github("fasterius/seqCAT")
```

You may also install `seqCAT` using [Conda][13]:

```bash
conda install -c bioconda bioconductor-seqcat
```

To list the versions of `seqCAT` available on Conda, you can use the `search`
functionality:

```bash
conda search -c bioconda bioconductor-seqcat
```

## Usage

The general workflow of `seqCAT` consists of three steps:

    1.  Creation of SNV profiles
    2.  Comparisons of SNV profiles
    3.  Evaluation of profile comparisons

```r
# Load the package
library("seqCAT")

# Path to the example VCF file
vcf <- system.file("extdata", "example.vcf.gz", package = "seqCAT")

# Create SNV profiles
hct116 <- create_profile(vcf, "HCT116")
hke3 <- create_profile(vcf, "HKE3")
rko <- create_profile(vcf, "RKO")

# Compare all profiles to each other
profiles <- list(hct116, hke3, rko)
comparisons <- compare_many(profiles)

# Create an heatmap of comparisons and their similarity scores
plot_heatmap(comparisons[[1]])
```
<p align="center">
    <img src="man/figures/README_example_1.png", alt="Example heatmap"/>
</p>

For more detailed instructions on how to use `seqCAT`, please see the
[vignette][14].

## Citation

If you are using `seqCAT` to analyse your data, please cite the
following article:

> **seqCAT: a Bioconductor R-package for variant analysis of high throughput**
> **sequencing data**
> <br/> Fasterius E. and Al-Khalili Szigyarto C.
> <br/> *F1000Research* (2018), 7:1466
> <br/> https://f1000research.com/articles/7-1466

## License

The `seqCAT` package is released with a MIT licence and is a free software: you
may redistribute and/or modify it under the terms of the license. For more
information, please see the `LICENCE` file.

[1]: https://anaconda.org/bioconda/bioconductor-seqcat/badges/version.svg
[2]: https://anaconda.org/bioconda/bioconductor-seqcat
[3]: https://img.shields.io/badge/License-MIT-blue.svg
[4]: https://opensource.org/licenses/MIT
[5]: https://travis-ci.org/fasterius/seqCAT.svg?branch=master
[6]: https://travis-ci.org/fasterius/seqCAT
[7]: https://coveralls.io/repos/github/fasterius/seqCAT/badge.svg?branch=master
[8]: https://coveralls.io/github/fasterius/seqCAT?branch=master

[9]: http://www.internationalgenome.org/wiki/Analysis/variant-call-format
[10]: https://software.broadinstitute.org/gatk/
[11]: http://snpeff.sourceforge.net/
[12]: https://bioconductor.org/packages/release/bioc/html/seqCAT.html
[13]: https://conda.io/en/latest/
[14]: https://bioconductor.org/packages/release/bioc/vignettes/seqCAT/inst/doc/seqCAT.html
