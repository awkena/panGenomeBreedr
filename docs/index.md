# panGenomeBreedr

`panGenomeBreedr` (`panGB`) is conceptualized to be a unified, crop
agnostic platform for pangenome-enabled breeding that follows
standardized conventions for natural or casual variant analysis using
pangenomes, marker design, and marker QC hypothesis testing (Figure 1).
It seeks to simplify and enhance the use of pangenome resources in
cultivar development.

![panGenomeBreedr marker design
workflow](reference/figures/workflow.png)

*Fig. 1.* Conceptual workflow of the panGenomeBreedr (panGB) package for
pangenome-enabled marker development. Using snpEff-annotated VCF files
organized within a SQLite database, panGB enables querying for variants
within candidate genes or user-defined genomic regions. It retrieves
mutation annotations and predicted impacts from snpEff to identify
putative causal variants (PCVs), which serve as the basis for designing
functional trait-predictive markers. The package supports
hypothesis-driven validation of these markers and also facilitates the
design of additional marker types, including precision-introgression and
background markers.

In its current development version, panGB provides customizable R
functions for **variant discovery from snpEff-annotated VCF files, KASP
marker design, and marker validation** (Steps 1–3 in Fig. 1).

To expand accessibility, `panGB` will include a user-friendly Shiny
application, allowing non-R users to leverage its core features without
requiring R programming experience.

The SNP Viewer tool by LGC Genomics is limited to Windows platforms and
lacks standardized conventions for visualizing positive controls in
marker validation, making it difficult for users to conclusively assess
marker performance. In contrast, `panGB` offers platform-independent
tools for hypothesis testing, quality control (QC), and validation of
KASP markers, addressing a key gap in existing visualization and
validation workflows.

Submit bug reports and feature suggestions, or track changes on the
[issues page](https://github.com/awkena/panGenomeBreedr/issues).

# Requirements

To use this package locally on a machine, the following R packages are
required:

- [ggplot2](https://CRAN.R-project.org/package=ggplot2)

- [gridExtra](https://CRAN.R-project.org/package=gridExtra)

- [reshape2](https://CRAN.R-project.org/package=reshape2)

- [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html)

- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)

- [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)

- [msa](https://bioconductor.org/packages/release/bioc/html/msa.html)

# Recommended packages

- [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html):
  Needed for package development and installation from GitHub on Windows
  PCs.

- [UpSetR](https://CRAN.R-project.org/package=UpSetR): Required for
  generating UpSet plots.

## Installation

First, ensure all existing packages are up to date.

You can install the development version of `panGenomeBreedr` from
[GitHub](https://github.com/awkena/panGenomeBreedr) with:

``` r
# Install panGenomeBreedr
if (!require("devtools")) install.packages("devtools")

devtools::install_github("awkena/panGenomeBreedr")
```

### Installing Bioconductor dependency packages

`panGB` depends on a list of Bioconductor packages that may not be
installed automatically alongside `panGB`. To manually install these
packages, use the code snippet below:

``` r
# Install and load required Bioconductor packages if not already installed
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Define required Bioconductor packages
bioc_packages <- c("rtracklayer", "Rsamtools", "msa", "IRanges",
                   "GenomicRanges", "BSgenome", "Biostrings")

# Install any missing Bioconductor packages
for (pkg in bioc_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    BiocManager::install(pkg, ask = FALSE, update = FALSE)

  }
}
```

# Current Functionality of `panGB`

`panGB` currently provides functionality for the following key tasks:

1.  **Variant discovery**  
    Identify variants within candidate genes or any user-defined genomic
    interval using *snpEff*-annotated VCF files.

2.  **KASP marker design**  
    Generate allele-specific markers targeting either causal variants or
    any variant of interest.

3.  **Marker validation and QC visualization**  
    Produce quality control plots and perform hypothesis-driven
    evaluations to assess marker reliability.

4.  **Decision-support for trait introgression**  
    Guide marker-assisted backcrossing by profiling foreground,
    background, and precision-introgression markers to support selection
    decisions.

👉 For a full tutorial and worked example, check out the
[panGenomeBreedr Workflow
vignette](https://awkena.github.io/panGenomeBreedr/articles/panGenomeBreedr_Workflows.html).

## Troubleshooting

If the package does not run as expected, check the following:

- Was the package properly installed?

- Do you have the required dependencies installed?

- Were any warnings or error messages returned during package
  installation?

- Are all packages up to date before installing panGB?

# Authors and contributors

- [Alexander Wireko Kena](https://github.com/awkena)

- [Israel Tawiah Tetteh](https://github.com/Israel-Tetteh)

- [Cruet Burgos](https://www.morrislab.org/people/clara-cruet-burgos)

- [Fanna Maina](https://www.morrislab.org/people/fanna-maina)

- [Linly Banda](https://www.biofortificationlab.org/people/linly-banda)

- [Jacques
  Faye](https://sites.google.com/site/morrislaboratory/people/jacques-faye)

- [Benjamin
  Annor](https://webapps.knust.edu.gh/staff/dirsearch/profile/summary/093520aa7216.html)

- [Terry
  Felderhoff](https://www.agronomy.k-state.edu/about/people/faculty/felderhoff-terry/)

- [Geoffrey Preston
  Morris](https://www.morrislab.org/people/geoff-morris)

# License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

# Support and Feedback

For support and submission of feedback, email the maintainer **Alexander
Kena, PhD** at <alex.kena24@gmail.com>
