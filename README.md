
<!-- README.md is generated from README.Rmd. Please edit that file -->

# panGenomeBreedr <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
<!-- badges: end -->

`panGenomeBreedr` (`panGB`) is conceptualized to be a unified, crop
agnostic platform for pangenome-enabled breeding that follows
standardized conventions for natural or casual variant analysis using
pangenomes, marker design, and marker QC hypothesis testing (Figure 1).
It seeks to simplify and enhance the use of pangenome resources in
cultivar development.

| <img src='man/figures/workflow.png' align="center" style="width: 800px;" /> |
|:--:|
| *Fig. 1. Conceptual workflow of the panGenomeBreedr (panGB) package for pangenome-enabled marker development. Using snpEff-annotated VCF files organized within a SQLite database, panGB enables querying for variants within candidate genes or user-defined genomic regions. It retrieves mutation annotations and predicted impacts from snpEff to identify putative causal variants (PCVs), which serve as the basis for designing functional trait-predictive markers. The package supports hypothesis-driven validation of these markers and also facilitates the design of additional marker types, including precision-introgression and background markers.* |

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

# Table of contents

- [Requirements](#requirements)
- [Recommended packages](#recommended-packages)  
- [Installation](#installation)
- [Current Functionality of panGB](#current-functionality-of-pangb)
- [Usage](#usage)
- [Examples](#examples)
  - [Variant Discovery](#variant-discovery)
    - [Pangenome Data and Database
      Rationale](#pangenome-data-and-database-rationale)
    - [Recommended Schema for the SQLite
      Database](#recommended-schema-for-the-sqlite-database)
    - [Database Creation](#database-creation)
  - [KASP Marker Design](#kasp-marker-design)
  - [KASP Marker Validation](#kasp-marker-validation)
- [Other Breeder-Centered Functionalities in
  panGB](#other-breeder-centered-functionalities-in-pangb)
  - [Creating Heatmaps with `panGB`](#creating-heatmaps-with-pangb)
  - [Trait Introgression Hypothesis
    Testing](#trait-introgression-hypothesis-testing)
  - [Decision Support for Marker-Assisted Backcrossing in
    `panGB`](#decision-support-for-marker-assisted-backcrossing-in-pangb)
  - [Weighted RPP computation in
    panGB](#weighted-rpp-computation-in-pangb)
- [Troubleshooting](#troubleshooting)
- [Authors and contributors](#authors-and-contributors)
- [License](#license)
- [Support and Feedback](##support-and-feedback)

## Requirements

To run this package locally on a machine, the following R packages are
required:

- [ggplot2](https://ggplot2.tidyverse.org): Elegant Graphics for Data
  Analysis.

- [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html):
  Miscellaneous Functions for “Grid” Graphics.

- [utils](https://www.rdocumentation.org/packages/utils/versions/3.6.2):
  The R Utils Package.

- [stats](https://www.rdocumentation.org/packages/stats/versions/3.6.2)

- [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

- [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)

- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)

- [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)

- [msa](https://bioconductor.org/packages/release/bioc/html/msa.html)

## Recommended packages

- [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html):
  Needed for package development and installation from GitHub on Windows
  PCs.

- [rmarkdown](https://CRAN.R-project.org/package=rmarkdown): When
  installed, display of the project’s README.md will be rendered with R
  Markdown.

## Installation

First, ensure all existing packages are up to date.

You can install the development version of `panGenomeBreedr` from
[GitHub](https://github.com/awkena/panGenomeBreedr) with:

``` r
# Install panGenomeBreedr
if (!require("devtools")) install.packages("devtools")

devtools::install_github("awkena/panGenomeBreedr", upgrade = TRUE)
```

### Installing Bioconductor dependency packages

`panGB` depends on a list of Bioconductor packages that may not be
installed automatically alongside `panGB`. To manually install these
packages, use the code snippet below:

``` r
# Install and load required Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

  BiocManager::install(c("Bsgenome,
                         "Biostrings",
                         "GenomicRanges",
                         "IRanges",
                         "msa"))
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

# Usage

## Examples

Here, we provide examples on how to use `panGB` to design a KASP marker
based on a causal variant, as well as marker validation for any KASP
marker.

## Variant Discovery

### Pangenome Data and Database Rationale

The examples used in this documentation are based on **sorghum pangenome
resources** derived from whole-genome resequencing data of **1,676
sorghum lines**. Variant calling was performed using version **v5.1** of
the **BTx623** reference genome. The resulting **SNP** and **INDEL**
variants were functionally annotated using **snpEff**.

Direct querying of **snpEff-annotated VCF files** from R is often
**computationally slow and inefficient**, especially with large
datasets. To overcome this limitation, we built a **SQLite database**
that stores the variants, annotations, and genotypes in normalized
tables. This structure allows for **fast and flexible access** to
relevant data, supporting workflows for trait-predictive marker
discovery

We **strongly recommend the creation of similar databases for other
crops**. The SQLite format offers a compact, portable, and queryable
representation of pangenome-derived variant data, significantly
improving performance and reproducibility in genomic analysis pipelines.

A compressed format of the SQLite database for sorghum can be downloaded
[here](https://drive.google.com/file/d/1L4S7_ZGeFyu_bA7rRsmTpf9V8__VLB-R/view?usp=sharing).

### Recommended Schema for the SQLite Database

The SQLite database contains the following three key tables:

#### `variants`

Stores core variant information extracted from the VCF.

| Column         | Description                        |
|----------------|------------------------------------|
| `variant_id`   | Unique variant identifier          |
| `chrom`        | Chromosome name                    |
| `pos`          | Genomic position (1-based)         |
| `ref`          | Reference allele                   |
| `alt`          | Alternate allele                   |
| `variant_type` | Type of variant (e.g., SNP, indel) |

#### `annotations`

Contains functional annotations from **snpEff**, typically including
predicted effects, gene names, and functional categories.

| Column         | Description                                |
|----------------|--------------------------------------------|
| `variant_id`   | Foreign key linking to `variants`          |
| `gene_name`    | Sorghum gene ID (e.g., “Sobic.005G213600”) |
| `effect`       | Type of effect (e.g., missense_variant)    |
| `impact`       | snpEff predicted impact (e.g., HIGH)       |
| `feature_type` | Type of annotated feature (e.g., exon)     |
| …              | Additional snpEff annotation fields        |

#### `genotypes`

Stores genotype calls per sample for each variant.

| Column       | Description                         |
|--------------|-------------------------------------|
| `variant_id` | Foreign key linking to `variants`   |
| `chrom`      | Chromosome name                     |
| `pos`        | Genomic position                    |
| `sample1`    | Genotype for Sample 1 (e.g., “A:G”) |
| `sample2`    | Genotype for Sample 2 (e.g., “A:A”) |
| …            | Genotypes for other samples         |

### Database Creation

We generated the SQLite database a custom workflow that 1. Parses a
multi-sample VCF file annotated by snpEff, 2. Extracts variant,
annotation, and genotype data, 3. Writes the data into normalized
relational tables (`variants`, `annotations`, `genotypes`).

A prebuilt mini example database (`mini_sorghum_variant_vcf.db.gz`) is
included in the `extdata/` folder of the package.

### Query Variant Tables from SQLite Database

The `query_db()` function allows users to query specific tables within a
panGenomeBreedr-formatted SQLite database for variants, annotations, or
genotypes based on chromosome coordinates or candidate gene IDs.

This function retrieves records from one of the following tables in the
database:

- `variants`: Basic variant information (chromosome, position, REF/ALT
  alleles, etc.)
- `annotations`: Variant effect predictions (e.g., from snpEff)
- `genotypes`: Genotypic data across lines/samples plus the metadata of
  the variants.

Users can specify genomic coordinates (`chrom`, `start`, `end`) or a
candidate gene name (`gene_name`) to extract relevant entries.

If used correctly, the `query_db()` function returns a data frame
containing the filtered records from the selected table.

``` r
library(panGenomeBreedr)

# This example uses the mini SQLite database included in the package.
# Define a temporary path and decompress the example database
path <- tempdir()
mini_db <- system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
                       package = "panGenomeBreedr", mustWork = TRUE)
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)

# Query VCF genotypes within the genomic range: Chr05:75104537-75106403
gt_region <- query_db(db_path = mini_db_path,
                      chrom = "Chr05",
                      start = 75104537,
                      end = 75106403,
                      table_name = "genotypes")

# Query snpEff annotations within a candidate gene
annota_region <- query_db(db_path = mini_db_path,
                          chrom = "Chr05",
                          start = 75104537,
                          end = 75106403,
                          table_name = "annotations",
                          gene_name = "Sobic.005G213600")

# Clean up temporary files
unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE), 
       recursive = TRUE)

knitr::kable(gt_region[1:5, 1:10], 
             caption = 'Table 1: Queried genotypes for varaints from the SQLite database.', 
             format = 'html', 
             padding = 0, 
             booktabs = TRUE, 
             escape = FALSE)
```

<table>
<caption>
Table 1: Queried genotypes for varaints from the SQLite database.
</caption>
<thead>
<tr>
<th style="text-align:left;">
variant_id
</th>
<th style="text-align:left;">
chrom
</th>
<th style="text-align:left;">
pos
</th>
<th style="text-align:left;">
variant_type
</th>
<th style="text-align:left;">
ref
</th>
<th style="text-align:left;">
alt
</th>
<th style="text-align:left;">
IDMM
</th>
<th style="text-align:left;">
ISGC
</th>
<th style="text-align:left;">
ISGK
</th>
<th style="text-align:left;">
ISHC
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
INDEL_Chr05_75104541
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104541
</td>
<td style="text-align:left;">
INDEL
</td>
<td style="text-align:left;">
T
</td>
<td style="text-align:left;">
TGAC
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
</tr>
<tr>
<td style="text-align:left;">
SNP_Chr05_75104557
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104557
</td>
<td style="text-align:left;">
SNP
</td>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
T
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
</tr>
<tr>
<td style="text-align:left;">
SNP_Chr05_75104560
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104560
</td>
<td style="text-align:left;">
SNP
</td>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
T
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
</tr>
<tr>
<td style="text-align:left;">
INDEL_Chr05_75104564
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104564
</td>
<td style="text-align:left;">
INDEL
</td>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
CA
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
</tr>
<tr>
<td style="text-align:left;">
SNP_Chr05_75104568
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104568
</td>
<td style="text-align:left;">
SNP
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
T
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
<td style="text-align:left;">
0\|0
</td>
</tr>
</tbody>
</table>

``` r

knitr::kable(annota_region[1:5,], 
             caption = 'Table 2: Queried annotations for variants from the SQLite database.', 
             format = 'html', 
             padding = 0, 
             booktabs = TRUE, 
             escape = FALSE)
```

<table>
<caption>
Table 2: Queried annotations for variants from the SQLite database.
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
variant_id
</th>
<th style="text-align:left;">
allele
</th>
<th style="text-align:left;">
annotation
</th>
<th style="text-align:left;">
impact
</th>
<th style="text-align:left;">
gene_name
</th>
<th style="text-align:left;">
gene_id
</th>
<th style="text-align:left;">
feature_type
</th>
<th style="text-align:left;">
feature_id
</th>
<th style="text-align:left;">
transcript_biotype
</th>
<th style="text-align:left;">
rank
</th>
<th style="text-align:left;">
HGVS_c
</th>
<th style="text-align:left;">
HGVS_p
</th>
<th style="text-align:left;">
chrom
</th>
<th style="text-align:right;">
pos
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
INDEL_Chr05_75104541
</td>
<td style="text-align:left;">
TGAC
</td>
<td style="text-align:left;">
3_prime_UTR_variant
</td>
<td style="text-align:left;">
MODIFIER
</td>
<td style="text-align:left;">
Sobic.005G213600
</td>
<td style="text-align:left;">
Sobic.005G213600.v5.1
</td>
<td style="text-align:left;">
transcript
</td>
<td style="text-align:left;">
Sobic.005G213600.1.v5.1
</td>
<td style="text-align:left;">
protein_coding
</td>
<td style="text-align:left;">
2/2
</td>
<td style="text-align:left;">
c.*322\_*324dupGTC
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:right;">
75104541
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
SNP_Chr05_75104557
</td>
<td style="text-align:left;">
T
</td>
<td style="text-align:left;">
3_prime_UTR_variant
</td>
<td style="text-align:left;">
MODIFIER
</td>
<td style="text-align:left;">
Sobic.005G213600
</td>
<td style="text-align:left;">
Sobic.005G213600.v5.1
</td>
<td style="text-align:left;">
transcript
</td>
<td style="text-align:left;">
Sobic.005G213600.1.v5.1
</td>
<td style="text-align:left;">
protein_coding
</td>
<td style="text-align:left;">
2/2
</td>
<td style="text-align:left;">
c.\*309G\>A
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:right;">
75104557
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
SNP_Chr05_75104560
</td>
<td style="text-align:left;">
T
</td>
<td style="text-align:left;">
3_prime_UTR_variant
</td>
<td style="text-align:left;">
MODIFIER
</td>
<td style="text-align:left;">
Sobic.005G213600
</td>
<td style="text-align:left;">
Sobic.005G213600.v5.1
</td>
<td style="text-align:left;">
transcript
</td>
<td style="text-align:left;">
Sobic.005G213600.1.v5.1
</td>
<td style="text-align:left;">
protein_coding
</td>
<td style="text-align:left;">
2/2
</td>
<td style="text-align:left;">
c.\*306G\>A
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:right;">
75104560
</td>
</tr>
<tr>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
INDEL_Chr05_75104564
</td>
<td style="text-align:left;">
CA
</td>
<td style="text-align:left;">
3_prime_UTR_variant
</td>
<td style="text-align:left;">
MODIFIER
</td>
<td style="text-align:left;">
Sobic.005G213600
</td>
<td style="text-align:left;">
Sobic.005G213600.v5.1
</td>
<td style="text-align:left;">
transcript
</td>
<td style="text-align:left;">
Sobic.005G213600.1.v5.1
</td>
<td style="text-align:left;">
protein_coding
</td>
<td style="text-align:left;">
2/2
</td>
<td style="text-align:left;">
c.\*301dupT
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:right;">
75104564
</td>
</tr>
<tr>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
SNP_Chr05_75104568
</td>
<td style="text-align:left;">
T
</td>
<td style="text-align:left;">
3_prime_UTR_variant
</td>
<td style="text-align:left;">
MODIFIER
</td>
<td style="text-align:left;">
Sobic.005G213600
</td>
<td style="text-align:left;">
Sobic.005G213600.v5.1
</td>
<td style="text-align:left;">
transcript
</td>
<td style="text-align:left;">
Sobic.005G213600.1.v5.1
</td>
<td style="text-align:left;">
protein_coding
</td>
<td style="text-align:left;">
2/2
</td>
<td style="text-align:left;">
c.\*298C\>A
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:right;">
75104568
</td>
</tr>
</tbody>
</table>

The `query_db()` function has the following key arguments: - `db_path`:
*(character)* Path to the SQLite database file. - `table_name`:
*(character)* Name of the table to query (`"variants"`, `"annotations"`,
or `"genotypes"`). - `chrom`: *(character)* Chromosome name to filter by
(e.g., `"Chr05"`). - `start`, `end`: *(numeric)* Start and end genomic
coordinates for the query region. - `gene_name`: *(character)* Optional
Sobic gene ID (e.g., `"Sobic.005G213600"`) — used when querying the
`annotations` table to filter entries by gene.

### Filter Variants by Allele Frequency in a Genomic Region

The `query_by_af()` function allows users to extract variants from a
SQLite database based on **alternate allele frequency thresholds**
within a specified genomic region.

This is particularly useful for identifying **polymorphic sites** within
candidate gene regions or windows of interest that meet desired minor
allele frequency (MAF) thresholds for marker development.

An example usage for the `query_by_af()` function is shown in the code
snippet below:

``` r
library(panGenomeBreedr)

# Define temporary directory and decompress demo database
path <- tempdir()
mini_db <- system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
                       package = "panGenomeBreedr", mustWork = TRUE)
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)

# Filter variants with alt allele frequency between 1% and 99% in a defined region
filter_af <- query_by_af(db_path = mini_db_path,
                         min_af = 0.01,
                         max_af = 0.99,
                         chrom = "Chr05",
                         start = 75104537,
                         end = 75106403)

head(filter_af)
#>              variant_id chrom      pos    ref_af     alt_af
#> 2    SNP_Chr05_75104557 Chr05 75104557 0.8905131 0.10948687
#> 3    SNP_Chr05_75104560 Chr05 75104560 0.8896181 0.11038186
#> 4  INDEL_Chr05_75104564 Chr05 75104564 0.8854415 0.11455847
#> 7  INDEL_Chr05_75104573 Chr05 75104573 0.8848449 0.11515513
#> 9  INDEL_Chr05_75104585 Chr05 75104585 0.8836516 0.11634845
#> 10   SNP_Chr05_75104591 Chr05 75104591 0.9642005 0.03579952

# Clean up temporary files
unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE), 
       recursive = TRUE)

knitr::kable(filter_af[1:5,], 
             caption = 'Table 3: Filtered variants from the SQLite database.', 
             format = 'html', 
             padding = 0, 
             booktabs = TRUE, 
             escape = FALSE)
```

<table>
<caption>
Table 3: Filtered variants from the SQLite database.
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
variant_id
</th>
<th style="text-align:left;">
chrom
</th>
<th style="text-align:left;">
pos
</th>
<th style="text-align:right;">
ref_af
</th>
<th style="text-align:right;">
alt_af
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
SNP_Chr05_75104557
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104557
</td>
<td style="text-align:right;">
0.8905131
</td>
<td style="text-align:right;">
0.1094869
</td>
</tr>
<tr>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
SNP_Chr05_75104560
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104560
</td>
<td style="text-align:right;">
0.8896181
</td>
<td style="text-align:right;">
0.1103819
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
INDEL_Chr05_75104564
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104564
</td>
<td style="text-align:right;">
0.8854415
</td>
<td style="text-align:right;">
0.1145585
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
INDEL_Chr05_75104573
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104573
</td>
<td style="text-align:right;">
0.8848449
</td>
<td style="text-align:right;">
0.1151551
</td>
</tr>
<tr>
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
INDEL_Chr05_75104585
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:left;">
75104585
</td>
<td style="text-align:right;">
0.8836516
</td>
<td style="text-align:right;">
0.1163484
</td>
</tr>
</tbody>
</table>

The `query_by_af()` function has the following input parameters:

| Argument  | Type        | Description                               |
|-----------|-------------|-------------------------------------------|
| `db_path` | `character` | Path to the SQLite database (`.db` file). |
| `min_af`  | `numeric`   | Minimum allele frequency (default: `0`).  |
| `max_af`  | `numeric`   | Maximum allele frequency (default: `1`).  |
| `chrom`   | `character` | Chromosome name (e.g., `"Chr05"`).        |
| `start`   | `numeric`   | Start coordinate of the region.           |
| `end`     | `numeric`   | End coordinate of the region.             |

### Summarize SnpEff Annotation and Impact Distributions by Variant Type

The `query_ann_summary()` function provides a convenient way to
summarize the distribution of **SnpEff annotations** and **impact
categories** across variant types (e.g., SNPs, indels) within a defined
genomic region.

This function enables users to quickly assess the types and functional
implications of variants located within candidate genes or genomic
intervals of interest.

``` r
library(panGenomeBreedr)

# Prepare test database
path <- tempdir()
mini_db <- system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
                       package = "panGenomeBreedr", mustWork = TRUE)
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)

# Run annotation summary for region Chr05:75104537-75106403
ann_summary <- query_ann_summary(db_path = mini_db_path,
                                 chrom = "Chr05",
                                 start = 75104537,
                                 end = 75106403)

# View summaries
head(ann_summary$annotation_summary)
#>                 annotation variant_type count
#> 24   upstream_gene_variant          SNP   149
#> 12   upstream_gene_variant        INDEL    53
#> 19 downstream_gene_variant          SNP    46
#> 22        missense_variant          SNP    29
#> 23      synonymous_variant          SNP    21
#> 13     3_prime_UTR_variant          SNP    18
head(ann_summary$impact_summary)
#>     impact variant_type count
#> 8 MODIFIER          SNP   220
#> 4 MODIFIER        INDEL    82
#> 7 MODERATE          SNP    29
#> 6      LOW          SNP    21
#> 1     HIGH        INDEL     6
#> 3 MODERATE        INDEL     5

# Clean up
unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE), 
       recursive = TRUE)
```

The `query_ann_summary()` function has the following input parameters:

| Argument | Type | Description |
|----|----|----|
| `db_path` | `character` | Path to the SQLite database. |
| `chrom` | `character` | Chromosome name (e.g., `"Chr05"`). |
| `start`, `end` | `numeric` | Start and end coordinates of the target genomic region. |
| `annotations_table` | `character` | Name of the table storing snpEff annotations (default: `"annotations"`). |
| `variants_table` | `character` | Name of the table storing variant data (default: `"variants"`). |

The `query_ann_summary()` function returns a `list` with the following
elements: - `annotation_summary`: Data frame summarizing the count of
each SnpEff annotation grouped by variant type. - `impact_summary`: Data
frame summarizing the count of each SnpEff impact level (e.g., HIGH,
MODERATE) grouped by variant type. - `variant_type_totals`: Total count
of variants in the region grouped by variant type.

**The annotation summary shows that there are six (6) INDEL variants
with a HIGH impact on protein function.** To see these variants, we need
to use the `query_by_impact()` function, as shown below:

``` r

# Prepare test database
path <- tempdir()
mini_db <- system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
                       package = "panGenomeBreedr", mustWork = TRUE)
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)

# Extract low impact variant for a region or gene
high_variants <- query_by_impact(db_path = mini_db_path,
                                impact_level = 'high',
                                chrom = "Chr05",
                                end = 75106403)

# View summaries
head(high_variants)
#>             variant_id chrom      pos   ref   alt qual filter variant_type
#> 1 INDEL_Chr05_75104881 Chr05 75104881     G GTCGA    .   PASS        INDEL
#> 2 INDEL_Chr05_75105587 Chr05 75105587    GC     G    .   PASS        INDEL
#> 3 INDEL_Chr05_75105598 Chr05 75105598    GT     G    .   PASS        INDEL
#> 4 INDEL_Chr05_75106156 Chr05 75106156 CGTAT     C    .   PASS        INDEL
#> 5 INDEL_Chr05_75106295 Chr05 75106295     A   ATC    .   PASS        INDEL
#> 6 INDEL_Chr05_75106325 Chr05 75106325     G   GTA    .   PASS        INDEL
#>   allele         annotation impact        gene_name               gene_id
#> 1  GTCGA frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 2      G frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 3      G frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 4      C frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 5    ATC frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#> 6    GTA frameshift_variant   HIGH Sobic.005G213600 Sobic.005G213600.v5.1
#>   feature_type              feature_id transcript_biotype rank
#> 1   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 2   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 3   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 4   transcript Sobic.005G213600.1.v5.1     protein_coding  2/2
#> 5   transcript Sobic.005G213600.1.v5.1     protein_coding  1/2
#> 6   transcript Sobic.005G213600.1.v5.1     protein_coding  1/2
#>               HGVS_c     HGVS_p
#> 1 c.1343_1344insTCGA p.Ser449fs
#> 2          c.637delG p.Ala213fs
#> 3          c.626delA p.Asp209fs
#> 4     c.65_68delATAC  p.His22fs
#> 5       c.38_39dupGA  p.Ser14fs
#> 6         c.8_9dupTA   p.Gln4fs

# Clean up
unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE), 
       recursive = TRUE)
```

## KASP Marker Design

The `kasp_marker_design()` function enables the design of **KASP
(Kompetitive Allele Specific PCR)** markers from identified putative
causal variants. It supports SNP, insertion, and deletion variants using
VCF genotype data and a reference genome to generate Intertek-compatible
marker information, including upstream and downstream polymorphic
context.

This function automates the extraction of flanking sequences and
polymorphic variants surrounding a focal variant and generates:

- Intertek-ready marker submission metadata
- DNA sequence alignment for visual inspection of marker context
- An optional publication-ready alignment plot in PDF format

The vcf file must contain the variant ID, Chromosome ID, Position, REF
and ALT alleles, as well as the genotype data for samples, as shown in
Table 1:

    #>             variant_id chrom      pos variant_type ref  alt IDMM ISGC ISGK ISHC
    #> 1 INDEL_Chr05_75104541 Chr05 75104541        INDEL   T TGAC  0|0  0|0  0|0  0|0
    #> 2   SNP_Chr05_75104557 Chr05 75104557          SNP   C    T  0|0  0|0  0|0  0|0
    #> 3   SNP_Chr05_75104560 Chr05 75104560          SNP   C    T  0|0  0|0  0|0  0|0
    #> 4 INDEL_Chr05_75104564 Chr05 75104564        INDEL   C   CA  0|0  0|0  0|0  0|0
    #> 5   SNP_Chr05_75104568 Chr05 75104568          SNP   G    T  0|0  0|0  0|0  0|0
    #> 6   SNP_Chr05_75104569 Chr05 75104569          SNP   C    A  0|0  0|0  0|0  0|0
    #>   ISHJ ISHK ISHS ISHT ISHU ISIC ISID ISIK ISIL ISIS ISIT ISIU ISJB ISJC ISJJ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   ISJK ISJS ISKB ISQC ISQE ISQF ISQG ISQH ISTM IUEB IUEC IUED IUEE IUEF IUEG
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  1|1
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  1|1
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IUEH IUEI IUEJ IUEK IUEL IUEM IUEN IUEP IUEQ IUER IUES IUET IUEU IUEW IUEX
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IUEY IUEZ IUFA IUFB IUFC IUFD IUFE IUFF IUFG IUFH IUFI IUFJ IUFK IUFL IUFN
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IUFS IUFT IUFU IUFW IUFX IUFY IUFZ IUGA IUGB IUGC IUGD IUGE IUGF IUGG IUGH
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IUGI IUGJ IUGK IUGL IUGM IUGN IUGP IUGQ IUGR IUGS IUGT IUGU IUGW IUGX IUGY
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  1|1  1|1
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  1|1  1|1
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  1|1  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IUGZ IUHA IUHB IUHC IUHD IUHE IUHF IUHG IUHH IUHI IUHJ IUHK IUHL IUHM IUHN
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|1  1|1  0|0  1|1  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|1  1|1  0|0  1|1  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|1  1|1  0|0  1|1  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IUHP IUHQ IUHR IUHS IUHT IUHU IUHW IXQW IXQX IXQY IXQZ IXRA IZKY IZLD IZLM
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZLN IZLP IZLQ IZLR IZLS IZLT IZLU IZLW IZLX IZLY IZLZ IZMA IZMB IZMC IZMD
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZME IZMF IZMG IZMH IZMI IZMJ IZMK IZML IZMM IZMN IZMP IZMQ IZMR IZMS IZMT
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZMU IZMX IZMY IZMZ IZNB IZND IZNE IZNF IZNG IZNH IZNI IZNJ IZNK IZNL IZNM
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  1|1  0|0
    #> 3  0|0  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  1|1  0|0
    #> 4  0|0  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  1|1  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZNN IZNP IZNQ IZNR IZNS IZNT IZNU IZNW IZNY IZNZ IZPA IZPB IZPC IZPD IZPE
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  1|1  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  1|1  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  1|1  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZPF IZPG IZPH IZPI IZPJ IZPK IZPL IZPM IZPN IZPP IZPQ IZPR IZPS IZPT IZPU
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 3  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 4  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZPW JBNU JBNW JBNX JBNY JBNZ JBPA JBPB JBPC JBPD JBPE JBPF JBPG JBPH JBPI
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  1|1
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  1|1
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   JBPJ JBPQ JBQM JBQU JBRM JBRZ JBSP IQRD JHIK IQRE ISJT IQRF JHIE JHIF JHGC
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   JHHY JHGK JHGJ JHGI JHGH JHGG IQRH IQRI IQRJ IQRK IQRL IQRM IQRN IQRP IQRQ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  1|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IQRR IQRS IQRT IQRU IQRW IQRX JHGF JHFQ IQRY IQRZ JHFP IQSA IQSB IQSC ISGD
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IQSD IQSE IQSF IQSG IQSH IQSI IQSJ IQSK IQSL JHID IQSY JHGL IQTB IQTE ISHD
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  1|1
    #> 3  0|0  0|0  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  1|1
    #> 4  0|0  0|0  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IQTF IQTG IQTH IQTI IQTJ IQTK IQTL ISHL JHGD JHIA IZGX IZGY JBQJ JBQK JBQL
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  1|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  1|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZHN JHGX IZGZ IZHA IZHB JBQN JBQP IZLG JBQQ IZHP JBQR JHGY JBQS IZHC JBQT
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZHD JBQW JBQX IZHE JBQY IZHF IZHG JBQZ JBRA JBRB JBRC IZHQ IZHR IZHS IZHT
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   JBRD IZHU IZHW JBRE IZHX IZHY JBRF JBRG JBRH JBRI JBRJ JBRK IZHH JBRL JBRN
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZHZ JBRP JBRQ IZIA IZIB IZIC IZID IZIE JBRR JBRS IZIF IZIG IZIH IZII IZIJ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZIK IZIL IZIM IZIN IZIP IZIQ IZIR IZIS IZIT JBRT IZIU IZIW IZIX IZIY IZIZ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZJA IZJB IZJC IZJD IZJE JBRU JBRW JBRX JBRY JBSA JBSB JBSC JBSD IZJF IZHI
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZHJ JBSE JBSF JBSG IZHK JBSH JBSI JBSJ IZJG JBSK JBSL JBSM IZJH IZJI IZJJ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   JBSN IZHL JBSQ IZJK IZJL IZHM JBSR IZJM IZJN IZJP JBSS JBST IZJQ IZJR IZJS
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZJT IZJU IZJW IZJX IZJY IZJZ IZKA INSR INSS INST INSU INSW INSX INSY INSZ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZKU INTB INTC INTD INTE INTF INTG INTH INTI INTJ IZKW INTL INTM INTN INTP
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INTQ INTR INTS INTT INTU INTW INTX INTY INTZ INUA INUB INUC INUD INUE INUF
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INUG IZLJ INUI INUJ INUK INUL INUM INUN INUP INUQ IZKX INUS INUT INUU INUW
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 3  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 4  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IQSU IQSW IQSX IQSZ IQTA IQTC IQTD IXMY IQSM IQSN IQSP IQSQ ISGL ISGU IQSR
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 3  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IQSS IQST JBPK JBPL JHFS JHEW JHFN INJP INJQ INJR INJS INJT IDYQ IFIT INUX
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   SRR8474782 SRR8474692 IFEE SAMEA12302664 INJU ISJD IFLR ISKI IFFZ INJW IDMN
    #> 1        0|0        0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2        0|0        0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3        0|0        0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4        0|0        0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5        0|0        0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6        0|0        0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INJX INJY INJZ INKA INKB INKC INKD ISKP INKE INKF IDYU IFMG ISJL IXPH IFGF
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   ISKW ISLC IFKN IDWG IFMN IFLC INKH ISLI IDWN IDWJ IDWI IFGA ISLP INKJ ISLW
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INKK IFFQ IFKS INKL IFKT INKM INKN IDTY INKP INKQ INKR INKS INKT INUY
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   SAMEA12302762 IFFJ IFGB IFFP INKU JHIC INKW INKX IDUX IFGW INUZ INKY INKZ
    #> 1           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0
    #> 4           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INLA INLB IFGR IFDU ISMI IFGY INWA INWB IFGK INLC INWC INWD INWE IDWT INWF
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1  0|0  0|0  1|1
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1  0|0  0|0  1|1
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1  0|0  0|0  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   ISMP ISMW JHGS INWH INWI IFLX INWJ ISJU INLE IDYG INWK INLF ISNB IDWK INLG
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IFML INLI INWL INWM ISIB INWN INWP INWQ INWR JHHU INWS INWT IDWD SRR8474623
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1        0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1        0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1        0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0
    #>   SRR8474811 SRR8474741 INWU ISKJ IECN SRR8474637 SRR8474714 IZPX SRR8474673
    #> 1        0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0        0|0
    #> 2        1|1        0|0  0|0  0|0  0|0        1|1        1|0  1|1        1|1
    #> 3        1|1        0|0  0|0  0|0  0|0        1|1        1|0  1|1        1|1
    #> 4        1|1        0|0  0|0  0|0  0|0        1|1        1|1  1|1        1|1
    #> 5        0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0        0|0
    #> 6        0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0        0|0
    #>   INLJ IZLC INWW INLK ISKX ISLD INWX IFDN IFDP ISGE IFGZ INWY IFJZ INLL INWZ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INXA IXQI INXC INXD INLM INXE INLN INLP IFJE IFJF INLQ INLR INLS IDUP IDZX
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  1|1  1|1  1|1  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  1|1  1|1  1|1  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  1|1  1|1  1|1  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INLT INXF INLU IFJG IFGI INLW INLX IFJH INLY ISLJ IFLM INLZ IXCU INMA INMB
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  0|0  0|0  0|0  0|0  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1  1|1
    #> 3  1|1  0|0  0|0  0|0  0|0  0|0  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1  1|1
    #> 4  1|1  0|0  0|0  0|0  0|0  0|0  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INMC INMD INME INMF ISLQ IFJJ IFJK IFJL IFJM IFJN INMG INMH IFHE IFEL INMI
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1  1|1  1|1  1|1
    #> 3  0|0  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1  1|1  1|1  1|1
    #> 4  0|0  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1  1|1  1|1  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INMJ INMK IFJP IFEM IFEN ISLX JHGT ISMD IFEP IFEQ IFER IFES IXPJ IFJQ IFJR
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  1|1  0|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1
    #> 3  1|1  1|1  0|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1
    #> 4  1|1  1|1  0|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  1|1  0|0  1|1
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   ISMJ IFJT SRR8474737 IXNM SRR8474732 IFET IFEU IFEW IFEX SRR8474731 INXH IFHF
    #> 1  0|0  0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0
    #> 2  0|0  1|1        0|0  0|0        0|0  0|0  1|1  0|0  0|0        1|1  0|0  0|0
    #> 3  0|0  1|1        0|0  0|0        0|0  0|0  1|1  0|0  0|0        1|1  0|0  0|0
    #> 4  0|0  1|1        0|0  0|0        0|0  0|0  1|1  0|0  0|0        1|1  0|0  0|0
    #> 5  0|0  0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0
    #> 6  0|0  0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0
    #>   IFEY IFEZ IFFA IFFB IFFC IFJX IFFD IFFE INXI IXML IFFF IFFG ISMQ IXPK INXJ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  1|1  1|1  1|1  1|1  1|1  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  1|1  1|1  1|1  1|1  1|1  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  1|1  1|1  1|1  1|1  1|1  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IFJY IFEH IZKP SRR8474634 SRR8474635 IXNN IFFH IXNP IFFI ISGM ISMX INML
    #> 1  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  1|1  0|0        0|0        0|0  1|1  1|1  1|1  1|1  0|0  0|0  0|0
    #> 3  1|1  1|1  0|0        0|0        0|0  1|1  1|1  1|1  1|1  0|0  0|0  0|0
    #> 4  1|1  1|1  0|0        0|0        0|0  1|1  1|1  1|1  1|1  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   SRR8474629 INMM INMN INMP INMQ IFEA IFEB IFEC INMR ISNC IFED IFLH IFLZ INXL
    #> 1        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2        0|0  0|0  0|0  0|0  0|0  0|0  0|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 3        0|0  0|0  0|0  0|0  0|0  0|0  0|1  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 4        0|0  0|0  0|0  0|0  0|0  0|0  1|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1
    #> 5        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   INXM IXQJ INXP INXQ IXQK IPBS IPBT SAMEA12302697 ISKK ISKR ISKY INMT ISLE
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  0|0  1|1  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  1|1  0|0  0|0  1|1  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  0|0  1|1  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0
    #>   ISLK ISLR JBPR ISLY IPBU IXRB IPBW IPBX IPBY IPBZ IPCA IQTM JHGM ISGW IDYP
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IXMM IXMN IXMP IPCC IXMQ IXMR IXMS IXMT IXMU IXMW IPCD IXMX IXPL IYYN IXNA
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IXPM IXNB IPCE IXPN IXPP IZKJ IXNC IXND IXPQ IXPR IXPS IXPT INMU IXPU IPCF
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZLA IPCG IXPW IXPX IXPY IXPZ INMW IPCH IXQA IXQB IZLH IXQC IPCI IPCJ INMX
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IFGN IDWL IDWM INMY IPCK IPCL ISHE SRR8474631 INMZ IFGG ISHM IPCM IPCN ISHW
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  1|1  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  1|1  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  1|1  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IPCP IFMB ISIE IXQL IPCR IZLI IPCT JHGA IPCW IXQM IDWU INNA IZKQ IPCZ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   SAMEA12302481 IFGX SAMEA12302509 SAMEA12302673 SAMEA12302666 SAMEA12302663
    #> 1           0|0  0|0           0|0           0|0           0|0           0|0
    #> 2           0|0  0|0           0|0           0|0           0|0           0|0
    #> 3           0|0  0|0           0|0           0|0           0|0           0|0
    #> 4           0|0  0|0           0|0           0|0           0|0           0|0
    #> 5           0|0  0|0           0|0           0|0           0|0           0|0
    #> 6           0|0  0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302643 SAMEA12302561 SAMEA12302632 SAMEA12302650 ISIM IZPZ
    #> 1           0|0           0|0           0|0           0|0  0|0  0|0
    #> 2           0|0           0|0           0|0           0|0  0|0  0|0
    #> 3           0|0           0|0           0|0           0|0  0|0  0|0
    #> 4           0|0           0|0           0|0           0|0  0|0  0|0
    #> 5           0|0           0|0           0|0           0|0  0|0  0|0
    #> 6           0|0           0|0           0|0           0|0  0|0  0|0
    #>   SAMEA12302629 SAMEA12302564 SAMEA12302562 SAMEA12302515 SAMEA12302655 IFGM
    #> 1           0|0           0|0           0|0           0|0           0|0  0|0
    #> 2           0|0           0|0           0|0           0|0           0|0  0|0
    #> 3           0|0           0|0           0|0           0|0           0|0  0|0
    #> 4           0|0           0|0           0|0           0|0           0|0  0|0
    #> 5           0|0           0|0           0|0           0|0           0|0  0|0
    #> 6           0|0           0|0           0|0           0|0           0|0  0|0
    #>   SAMEA12302599 SAMEA12302551 ISIW SAMEA12302544 SAMEA12302656 ISJE
    #> 1           0|0           0|0  0|0           0|0           0|0  0|0
    #> 2           0|0           0|0  0|0           0|0           1|1  0|0
    #> 3           0|0           0|0  0|0           0|0           1|1  0|0
    #> 4           0|0           0|0  0|0           0|0           1|1  0|0
    #> 5           0|0           0|0  0|0           0|0           0|0  0|0
    #> 6           0|0           0|0  0|0           0|0           0|0  0|0
    #>   SAMEA12302441 SAMEA12302474 SAMEA12302783 SAMEA12302630 SAMEA12302833
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           1|1
    #> 3           0|0           0|0           0|0           0|0           1|1
    #> 4           0|0           0|0           0|0           0|0           1|1
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302838 SAMEA12302644 SAMEA12302637 SAMEA12302819 SAMEA12302590
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302448 SAMEA12302794 SAMEA12302775 SAMEA12302486 SAMEA12302769
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302471 IFLW SAMEA12302594 IWNY SAMEA12302580 SAMEA12302674 ISJM
    #> 1           0|0  0|0           0|0  0|0           0|0           0|0  0|0
    #> 2           0|0  0|0           0|0  0|0           0|0           0|0  0|0
    #> 3           0|0  0|0           0|0  0|0           0|0           0|0  0|0
    #> 4           0|0  0|0           0|0  0|0           1|1           0|0  0|0
    #> 5           0|0  0|0           0|0  0|0           0|0           0|0  0|0
    #> 6           0|0  0|0           0|0  0|0           0|0           0|0  0|0
    #>   SAMEA12302577 SAMEA12302595 SAMEA12302659 INNB SAMEA12302772 SAMEA12302825
    #> 1           0|0           0|0           0|0  0|0           0|0           0|0
    #> 2           0|0           0|0           0|0  1|1           1|1           1|1
    #> 3           0|0           0|0           0|0  1|1           1|1           1|1
    #> 4           0|0           0|0           0|0  1|1           1|1           1|1
    #> 5           0|0           0|0           0|0  0|0           0|0           0|0
    #> 6           0|0           0|0           0|0  0|0           0|0           0|0
    #>   SAMEA12302649 SAMEA12302834 SAMEA12302651 SAMEA12302808 SAMEA12302791
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302497 SAMEA12302513 SAMEA12302798 SAMEA12302790 SAMEA12302727
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           1|1           1|1           1|1           0|0           0|0
    #> 3           1|1           1|1           1|1           0|0           0|0
    #> 4           1|1           1|1           1|1           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302712 SAMEA12302797 SAMEA12302741 SAMEA12302504 SAMEA12302698
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302611 ISJW SAMEA12302699 ISGF SAMEA12302793 ISGN SAMEA12302615
    #> 1           0|0  0|0           0|0  0|0           0|0  0|0           0|0
    #> 2           0|0  0|0           0|0  0|0           0|0  0|0           0|0
    #> 3           0|0  0|0           0|0  0|0           0|0  0|0           0|0
    #> 4           0|0  0|0           0|0  0|0           0|0  0|0           0|0
    #> 5           0|0  0|0           0|0  0|0           0|0  0|0           0|0
    #> 6           0|0  0|0           0|0  0|0           0|0  0|0           0|0
    #>   SAMEA12302636 SAMEA12302690 SAMEA12302689 SAMEA12302646 SAMEA12302635
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302618 SAMEA12302728 SAMEA12302638 SAMEA12302645 SAMEA12302613
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302639 SAMEA12302685 SAMEA12302809 SAMEA12302693 SAMEA12302778 INNC
    #> 1           0|0           0|0           0|0           0|0           0|0  0|0
    #> 2           0|0           0|0           1|1           0|0           0|0  0|0
    #> 3           0|0           0|0           1|1           0|0           0|0  0|0
    #> 4           0|0           0|0           1|1           0|0           0|0  0|0
    #> 5           0|0           0|0           0|0           0|0           0|0  0|0
    #> 6           0|0           0|0           0|0           0|0           0|0  0|0
    #>   SAMEA12302508 ISGX SAMEA12302462 SAMEA12302483 SAMEA12302597 SAMEA12302817
    #> 1           0|0  0|0           0|0           0|0           0|0           0|0
    #> 2           0|0  0|0           0|0           0|0           0|0           0|0
    #> 3           0|0  0|0           0|0           0|0           0|0           0|0
    #> 4           0|0  0|0           0|0           0|0           0|0           0|0
    #> 5           0|0  0|0           0|0           0|0           0|0           0|0
    #> 6           0|0  0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302660 ISHF SAMEA12302648 SAMEA12302748 SAMEA12302574 SAMEA12302572
    #> 1           0|0  0|0           0|0           0|0           0|0           0|0
    #> 2           1|1  0|0           0|0           0|0           0|0           0|0
    #> 3           1|1  0|0           0|0           0|0           0|0           0|0
    #> 4           1|1  1|1           0|0           0|0           0|0           0|0
    #> 5           0|0  0|0           0|0           0|0           0|0           0|0
    #> 6           0|0  0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302566 SAMEA12302686 ISHN SAMEA12302575 SAMEA12302576 SAMEA12302583
    #> 1           0|0           0|0  0|0           0|0           0|0           0|0
    #> 2           0|0           1|1  0|0           0|0           0|0           0|0
    #> 3           0|0           1|1  0|0           0|0           0|0           0|0
    #> 4           0|0           1|1  0|0           0|0           0|0           0|0
    #> 5           0|0           0|0  0|0           0|0           0|0           0|0
    #> 6           0|0           0|0  0|0           0|0           0|0           0|0
    #>   SAMEA12302739 SAMEA12302609 SAMEA12302687 SAMEA12302623 SAMEA12302617
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           1|1           0|0           0|0
    #> 3           0|0           0|0           1|1           0|0           0|0
    #> 4           0|0           0|0           1|1           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302642 SAMEA12302647 ISHX SAMEA12302803 SAMEA12302598 INND
    #> 1           0|0           0|0  0|0           0|0           0|0  0|0
    #> 2           0|0           0|0  0|0           1|1           1|1  1|1
    #> 3           0|0           0|0  0|0           1|1           1|1  1|1
    #> 4           0|0           0|0  0|0           1|1           1|1  1|1
    #> 5           0|0           0|0  0|0           0|0           0|0  0|0
    #> 6           0|0           0|0  0|0           0|0           0|0  0|0
    #>   SAMEA12302824 SAMEA12302604 SAMEA12302493 SAMEA12302477 ISIF SAMEA12302492
    #> 1           0|0           0|0           0|0           0|0  0|0           0|0
    #> 2           1|1           1|1           0|0           0|0  0|0           0|0
    #> 3           1|1           1|1           0|0           0|0  0|0           0|0
    #> 4           1|1           1|1           0|0           0|0  0|0           0|0
    #> 5           0|0           0|0           0|0           0|0  0|0           0|0
    #> 6           0|0           0|0           0|0           0|0  0|0           0|0
    #>   IZGS SAMEA12302457 SAMEA12302831 SAMEA12302801 SAMEA12302815 SAMEA12302567
    #> 1  0|0           0|0           0|0           0|0           0|0           0|0
    #> 2  0|0           0|0           0|0           0|0           0|0           0|0
    #> 3  0|0           0|0           0|0           0|0           0|0           0|0
    #> 4  0|0           0|0           0|0           0|0           0|0           0|0
    #> 5  0|0           0|0           0|0           0|0           0|0           0|0
    #> 6  0|0           0|0           0|0           0|0           0|0           0|0
    #>   ISIN ISIX ISJF IFHK ISJN INNE SAMEA12302500 IFMK IDYI IEBU IFMJ IDXC INNH
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  1|1  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  1|1  0|0  1|1  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  1|1  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IDYK IDWR SAMEA12302841 SAMEA12302675 IFGP SAMEA12302620 IDXE INNI
    #> 1  0|0  0|0           0|0           0|0  0|0           0|0  0|0  0|0
    #> 2  0|0  0|0           0|0           0|0  0|0           0|0  0|0  0|0
    #> 3  0|0  0|0           0|0           0|0  0|0           0|0  0|0  0|0
    #> 4  0|0  0|0           0|0           0|0  0|0           0|0  0|0  0|0
    #> 5  0|0  0|0           0|0           0|0  0|0           0|0  0|0  0|0
    #> 6  0|0  0|0           0|0           0|0  0|0           0|0  0|0  0|0
    #>   SAMEA12302750 SAMEA12302832 INNJ SAMEA12302746 ISME IDWQ SAMEA12302511 INNK
    #> 1           0|0           0|0  0|0           0|0  0|0  0|0           0|0  0|0
    #> 2           0|0           0|0  0|0           0|0  0|0  0|0           0|0  0|0
    #> 3           0|0           0|0  0|0           0|0  0|0  0|1           0|0  0|0
    #> 4           0|0           0|0  0|0           0|0  0|0  0|0           0|0  0|0
    #> 5           0|0           0|0  0|0           0|0  0|0  0|0           0|0  0|0
    #> 6           0|0           0|0  0|0           0|0  0|0  0|0           0|0  0|0
    #>   IPDA INNL IFLJ INNN ISMK IXNE IPDB ISMR ISMY ISND IFLK IFLL IDZN ISKL IDZP
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   SRR8474664 IDWF IDWP IEBN ISKS IDYM ISJX IQTN ISKZ IPDC IZKC IEDF IPDD IXQD
    #> 1        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IPDE IFEG IFEI IFEJ IZKD IFJD IFEK ISLF IZKE IZKK IZKF IPDF IXRC IZKG IUHX
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   JHGR SRR8474670 IPDH IQTP IQTQ IPDI IPDJ IQTR ISGG IWNU SAMEA12302560
    #> 1  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0
    #> 2  0|0        1|1  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0           0|0
    #> 3  0|0        1|1  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0           0|0
    #> 4  0|0        1|1  1|1  0|0  0|0  1|1  0|0  0|0  0|0  0|0           0|0
    #> 5  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0
    #> 6  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0
    #>   SAMEA12302553 ISGP IFDY IZNA IFGE JHGN IPDL IPDM IPDN JHGP IPDQ IPDR JHGB
    #> 1           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IFFK IFFL IPDS IPDT IFKA IFKB ISGY IFIX IFKP IFIY IFJC IFKY IFFU IFDQ IFFY
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IFDT IFFR IFDS IFGD IFDZ IFGL IFFT IFIZ IFFW IFGC IFKX IFFX IFKZ IFFS IXNQ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IZKZ IFFM IPDU IFKC IFKD IFKE IFKF IFKG IFKH ISLL IDUY IDMK IFLP ISLS ISLZ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 3  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   ISMF ISML IDWA ISMS ISMZ IDML ISNE ISKM IXNR IXNS IXNT IZKH IPDW ISKT IXNU
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IXNW IPDX IPDY IXNX IXNY IPDZ IXNZ IPEA IPEB IZKI IXPA IXPB IDWB IFLT IPEC
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IFMH IEAZ ISLA IFKI JHII IXPC IUHY IUHZ IUIA IXPD IXPE IXPF IXMF IXMG IXMH
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   ISHG IPED IPEE IPEF IPEG IPEH SRR8474581 IPEI IZKL IPEJ IXPG IPEK IXMI IPEL
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IPEM ISLG IPEN IPEP ISHP SAMEA12302725 SAMEA12302782 SAMEA12302717
    #> 1  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0
    #> 2  0|1  0|0  0|0  0|0  0|0           1|1           0|0           0|0
    #> 3  0|1  0|0  0|0  0|0  0|0           1|1           0|0           0|0
    #> 4  1|0  0|0  0|0  0|0  0|0           1|1           0|0           0|0
    #> 5  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0
    #> 6  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0
    #>   SAMEA12302701 ISHY SAMEA12302733 SAMEA12302811 SAMEA12302796 ISIG
    #> 1           0|0  0|0           0|0           0|0           0|0  0|0
    #> 2           0|0  1|1           0|0           0|0           1|0  0|0
    #> 3           0|0  1|1           0|0           0|0           1|0  0|0
    #> 4           0|0  1|1           0|0           0|0           0|1  0|0
    #> 5           0|0  0|0           0|0           0|0           0|0  0|0
    #> 6           0|0  0|0           0|0           0|0           0|0  0|0
    #>   SAMEA12302800 SAMEA12302702 SAMEA12302517 SAMEA12302723 IFHG SAMEA12302810
    #> 1           0|0           0|0           0|0           0|0  0|0           0|0
    #> 2           1|1           0|0           0|0           0|0  0|0           0|0
    #> 3           1|1           0|0           0|0           0|0  0|0           0|0
    #> 4           1|1           0|0           0|0           0|0  0|0           0|0
    #> 5           0|0           0|0           0|0           0|0  0|0           0|0
    #> 6           0|0           0|0           0|0           0|0  0|0           0|0
    #>   SAMEA12302764 SAMEA12302502 SAMEA12302473 SAMEA12302540 SAMEA12302476 ISIP
    #> 1           0|0           0|0           0|0           0|0           0|0  0|0
    #> 2           0|0           0|0           0|0           0|0           0|0  0|0
    #> 3           0|0           0|0           0|0           0|0           0|0  0|0
    #> 4           0|0           0|0           0|0           0|0           0|0  0|0
    #> 5           0|0           0|0           0|0           0|0           0|0  0|0
    #> 6           0|0           0|0           0|0           0|0           0|0  0|0
    #>   ISIY SAMEA12302628 SAMEA12302640 SAMEA12302489 SAMEA12302465 SAMEA12302732
    #> 1  0|0           0|0           0|0           0|0           0|0           0|0
    #> 2  0|0           0|0           0|0           0|0           0|0           0|0
    #> 3  0|0           0|0           0|0           0|0           0|0           0|0
    #> 4  0|0           0|0           0|0           0|0           0|0           0|0
    #> 5  0|0           0|0           0|0           0|0           0|0           0|0
    #> 6  0|0           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302641 ISJG SAMEA12302804 IFGJ ISJP SAMEA12302703 SAMEA12302472
    #> 1           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #> 2           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #> 3           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #> 4           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #> 5           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #> 6           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #>   SAMEA12302485 SAMEA12302546 SAMEA12302682 SAMEA12302654 SAMEA12302829 ISLM
    #> 1           0|0           0|0           0|0           0|0           0|0  0|0
    #> 2           0|0           0|0           1|1           0|0           0|0  0|0
    #> 3           0|0           0|0           1|1           0|0           0|0  0|0
    #> 4           0|0           0|0           1|1           0|0           0|0  0|0
    #> 5           0|0           0|0           0|0           0|0           0|0  0|0
    #> 6           0|0           0|0           0|0           0|0           0|0  0|0
    #>   SAMEA12302558 SAMEA12302827 JHHX IPEQ ISLT IZKB IPER ISMA ISMG IEAT ISMM ISMT
    #> 1           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IXNG IEAX IXQE IZKM ISNA IPES IXQF IXQG IZKN IPET IXQH JHFR ISJY JHHR JHGZ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  1|1  0|0  0|0  1|1  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  1|1  0|0  0|0  1|1  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  1|1  0|0  0|0  1|1  1|1  1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   JHHA JHHB JBPS IDWC JHHC IFLG IXNH JHHD JHFZ IZLK JHHE IXNI IXNJ IPEU IXNK
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IXNL JHHF JHHG IPEW IPEX JHHH IFLD JHHI JHHJ JHHK IFFN JHHL IFMP IPEY IPEZ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IDWW SAMEA12302466 SAMEA12302704 ISGH SAMEA12302619 SAMEA12302626
    #> 1  0|0           0|0           0|0  1|0           0|0           0|0
    #> 2  0|0           0|0           0|0  0|0           0|0           1|1
    #> 3  0|0           0|0           0|0  0|0           0|0           1|1
    #> 4  0|0           0|0           0|0  0|0           0|0           1|1
    #> 5  0|0           0|0           0|0  1|1           0|0           0|0
    #> 6  0|0           0|0           0|0  1|1           0|0           0|0
    #>   SAMEA12302463 SAMEA12302667 SAMEA12302802 IFMD SAMEA12302761 SAMEA12302781
    #> 1           0|0           0|0           0|0  0|0           0|0           0|0
    #> 2           0|0           0|0           0|0  0|0           0|0           0|0
    #> 3           0|0           0|0           0|0  0|0           0|0           0|0
    #> 4           0|0           0|0           0|0  0|0           0|0           0|0
    #> 5           0|0           0|0           0|0  0|0           0|0           0|0
    #> 6           0|0           0|0           0|0  0|0           0|0           0|0
    #>   SAMEA12302494 ISGQ SAMEA12302774 SAMEA12302785 SAMEA12302585 SAMEA12302816
    #> 1           0|0  0|0           0|0           0|0           0|0           0|0
    #> 2           0|0  1|1           0|0           0|0           0|0           0|0
    #> 3           0|0  1|1           0|0           0|0           0|0           0|0
    #> 4           0|0  1|1           0|0           0|0           0|0           0|0
    #> 5           0|0  0|0           0|0           0|0           0|0           0|0
    #> 6           0|0  0|0           0|0           0|0           0|0           0|0
    #>   IXNF SAMEA12302491 SAMEA12302765 SAMEA12302688 SAMEA12302548 SAMEA12302607
    #> 1  0|0           0|0           0|0           0|0           0|0           0|0
    #> 2  1|1           0|0           0|0           0|0           1|1           1|1
    #> 3  1|1           0|0           0|0           0|0           1|1           1|1
    #> 4  1|1           0|0           0|0           0|0           1|1           1|1
    #> 5  0|0           0|0           0|0           0|0           0|0           0|0
    #> 6  0|0           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302614 SAMEA12302537 SAMEA12302692 SAMEA12302658 SAMEA12302836
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           1|1           0|0           0|0
    #> 3           0|0           0|0           1|1           0|0           0|0
    #> 4           0|0           0|0           1|1           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302840 SAMEA12302569 IUIN SAMEA12302627 SAMEA12302496 SAMEA12302495
    #> 1           0|0           0|0  0|0           0|0           0|0           0|0
    #> 2           0|0           0|0  0|0           0|0           0|0           0|0
    #> 3           0|0           0|0  0|0           0|0           0|0           0|0
    #> 4           0|0           0|0  0|0           0|0           0|0           0|0
    #> 5           0|0           0|0  0|0           0|0           0|0           0|0
    #> 6           0|0           0|0  0|0           0|0           0|0           0|0
    #>   SAMEA12302522 IFHM SAMEA12302678 IPFA ISGZ JHHM ISHH IZKR JHHN ISHQ ISHZ ISIH
    #> 1           0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2           0|0  0|0           1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3           0|0  0|0           1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4           0|0  0|0           1|1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5           0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6           0|0  0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   SAMEA12302588 SRR8474587 IPFC IZLL IPFE IPFF IPFG IXQP JHIG IPFI IPFJ
    #> 1           0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2           0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 3           0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 4           0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0
    #> 5           0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6           0|0        0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   SAMEA12302839 IPFK IDXU IEBA IXQQ IPFM IPFN IPFP IPFQ IDWE JHHP IFKM IQTS
    #> 1           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6           0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IQTT JHHQ JBPT IQTU SAMEA12302830 SAMEA12302520 IFLY IFLE JBPU IFIW IFJA IFJB
    #> 1  0|0  0|0  0|0  0|0           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0           0|0           0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IFEF ISKN JBPW IFKQ IFKR ISKU IFLA SAMEA12302705 JHIJ SAMEA12302751 IFGU IFDW
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0
    #>   IFIU SAMEA12302498 IFGQ IFLB SAMEA12302552 ISIQ SAMEA12302612 JBPX IFKU IFKW
    #> 1  0|0           0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0  0|0
    #> 2  0|0           0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0  0|0
    #> 3  0|0           0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0  0|0
    #> 4  0|0           0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0  0|0
    #> 5  0|0           0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0  0|0
    #> 6  0|0           0|0  0|0  0|0           0|0  0|0           0|0  0|0  0|0  0|0
    #>   IQTW IZKS JHIP JHGU IQTY JHHT IQTZ JHHW IQUA IXQR IQUC IQUD IQUE IQUF JBPY
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0
    #>   JBPZ JBQA IQUG JBQB IFKJ JBQC IFKK IFKL JBQD JBQE IZKT JBQF IQUI IXQS IQUK
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   IQUL IQUM IXQT IZLE SRR8474602 JBQG SRR8474606 IQUQ IFDR SAMEA12302822 ISLB
    #> 1  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0  0|0           0|0  0|0
    #> 2  0|0  0|0  0|0  0|0        0|0  0|0        0|0  1|1  0|0           0|0  0|0
    #> 3  0|0  0|0  0|0  0|0        0|0  0|0        0|0  1|1  0|0           0|0  0|0
    #> 4  0|0  0|0  0|0  0|0        0|0  0|0        0|0  1|1  0|0           0|0  0|0
    #> 5  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0  0|0           0|0  0|0
    #> 6  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0  0|0           0|0  0|0
    #>   IFDX IQUR IFGT IFHI IFMR ISLN SAMEA12302694 SAMEA12302706 SAMEA12302449 JBQH
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0           0|0           0|0           0|0  0|0
    #>   ISIZ SAMEA12302707 SAMEA12302729 SAMEA12302718 SAMEA12302550 ISLU ISJH
    #> 1  0|0           0|0           0|0           0|0           0|0  0|0  0|0
    #> 2  0|0           0|0           0|0           0|0           0|0  0|0  0|0
    #> 3  0|0           0|0           0|0           0|0           0|0  0|0  0|0
    #> 4  0|0           0|0           0|0           0|0           0|0  0|0  0|0
    #> 5  0|0           0|0           0|0           0|0           0|0  0|0  0|0
    #> 6  0|0           0|0           0|0           0|0           0|0  0|0  0|0
    #>   SAMEA12302820 ISMB SAMEA12302749 SAMEA12302708 SAMEA12302758 SAMEA12302519
    #> 1           0|0  0|0           0|0           0|0           0|0           0|0
    #> 2           0|0  0|0           0|0           0|0           0|0           0|0
    #> 3           0|0  0|0           0|0           0|0           0|0           0|0
    #> 4           0|0  0|0           0|0           0|0           0|0           0|0
    #> 5           0|0  0|0           0|0           0|0           0|0           0|0
    #> 6           0|0  0|0           0|0           0|0           0|0           0|0
    #>   IFHB SAMEA12302726 IFHC SAMEA12302738 ISJQ SAMEA12302721 SAMEA12302709
    #> 1  0|0           0|0  0|0           0|0  0|0           0|0           0|0
    #> 2  0|0           0|0  0|0           0|0  0|0           0|0           0|0
    #> 3  0|0           0|0  0|0           0|0  0|0           0|0           0|0
    #> 4  0|0           0|0  0|0           0|0  0|0           0|0           0|0
    #> 5  0|0           0|0  0|0           0|0  0|0           0|0           0|0
    #> 6  0|0           0|0  0|0           0|0  0|0           0|0           0|0
    #>   SAMEA12302731 SAMEA12302760 SAMEA12302499 IFMC SAMEA12302521 SAMEA12302592
    #> 1           0|0           0|0           0|0  0|0           0|0           0|0
    #> 2           0|0           0|0           0|0  0|0           0|0           0|0
    #> 3           0|0           0|0           0|0  0|0           0|0           0|0
    #> 4           0|0           0|0           0|0  0|0           0|0           0|0
    #> 5           0|0           0|0           0|0  0|0           0|0           0|0
    #> 6           0|0           0|0           0|0  0|0           0|0           0|0
    #>   SAMEA12302805 SAMEA12302776 SAMEA12302823 JBQI SAMEA12302593 SAMEA12302737
    #> 1           0|0           0|0           0|0  0|0           0|0           0|0
    #> 2           0|0           0|0           0|0  0|0           0|0           0|0
    #> 3           0|0           0|0           0|0  0|0           0|0           0|0
    #> 4           0|0           0|0           0|0  0|0           0|0           0|0
    #> 5           0|0           0|0           0|0  0|0           0|0           0|0
    #> 6           0|0           0|0           0|0  0|0           0|0           0|0
    #>   SAMEA12302789 SAMEA12302754 SAMEA12302533 SAMEA12302557 ISGI SAMEA12302669
    #> 1           0|0           0|0           0|0           0|0  0|0           0|0
    #> 2           0|0           0|0           0|0           0|0  0|0           0|0
    #> 3           0|0           0|0           0|0           0|0  0|0           0|0
    #> 4           0|0           0|0           0|0           0|0  0|0           0|0
    #> 5           0|0           0|0           0|0           0|0  0|0           0|0
    #> 6           0|0           0|0           0|0           0|0  0|0           0|0
    #>   SAMEA12302442 SAMEA12302475 SAMEA12302505 ISGR SAMEA12302784 SAMEA12302443
    #> 1           0|0           0|0           0|0  0|0           0|0           0|0
    #> 2           0|0           0|0           0|0  0|0           0|0           0|0
    #> 3           0|0           0|0           0|0  0|0           0|0           0|0
    #> 4           0|0           0|0           0|0  0|0           0|0           0|0
    #> 5           0|0           0|0           0|0  0|0           0|0           0|0
    #> 6           0|0           0|0           0|0  0|0           0|0           0|0
    #>   SAMEA12302730 SAMEA12302745 IFHA SAMEA12302453 SAMEA12302719 ISHA IZPY
    #> 1           0|0           0|0  0|0           0|0           0|0  0|0  0|0
    #> 2           0|0           0|0  0|0           0|0           0|0  0|0  0|0
    #> 3           0|0           0|0  0|0           0|0           0|0  0|0  0|0
    #> 4           0|0           0|0  0|0           0|0           0|0  0|0  0|0
    #> 5           0|0           0|0  0|0           0|0           0|0  0|0  0|0
    #> 6           0|0           0|0  0|0           0|0           0|0  0|0  0|0
    #>   SAMEA12302523 ISHI SAMEA12302527 ISTN ISTP IWCK SAMEA12302633 SAMEA12302722
    #> 1           0|0  0|0           0|0  0|0  0|0  0|0           0|0           0|0
    #> 2           0|0  0|0           0|0  0|0  0|0  0|0           0|0           0|0
    #> 3           0|0  0|0           0|0  0|0  0|0  0|0           0|0           0|0
    #> 4           0|0  0|0           0|0  1|0  0|0  1|0           1|1           0|0
    #> 5           0|0  0|0           0|0  0|0  0|0  0|0           0|0           0|0
    #> 6           0|0  0|0           0|0  0|0  0|0  0|0           0|0           0|0
    #>   SAMEA12302488 IFLU ITIK SAMEA12302677 SAMEA12302755 SAMEA12302444
    #> 1           0|0  0|0  0|0           0|0           0|0           0|0
    #> 2           0|0  0|0  0|0           0|0           0|0           0|0
    #> 3           0|0  0|0  0|0           0|0           0|0           0|0
    #> 4           0|0  0|0  0|0           0|0           0|0           0|0
    #> 5           0|0  0|0  0|0           0|0           0|0           0|0
    #> 6           0|0  0|0  0|0           0|0           0|0           0|0
    #>   SAMEA12302532 ISIA SAMEA12302545 SAMEA12302724 SAMEA12302535 SAMEA12302602
    #> 1           0|0  0|0           0|0           0|0           0|0           0|0
    #> 2           0|0  1|1           0|0           0|0           0|0           0|0
    #> 3           0|0  1|1           0|0           0|0           0|0           0|0
    #> 4           0|0  1|1           0|0           0|0           0|0           0|0
    #> 5           0|0  0|0           0|0           0|0           0|0           0|0
    #> 6           0|0  0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302445 SAMEA12302742 SAMEA12302524 SAMEA12302525 SAMEA12302503
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           1|1
    #> 3           0|0           0|0           0|0           0|0           1|1
    #> 4           0|0           0|0           0|0           0|0           1|1
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302454 ISII SAMEA12302735 ITIL IFHH SAMEA12302766 SAMEA12302771
    #> 1           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #> 2           0|0  0|0           0|0  1|1  0|0           0|0           0|0
    #> 3           0|0  0|0           0|0  1|1  0|0           0|0           0|0
    #> 4           0|0  0|0           0|0  1|1  0|0           0|0           0|0
    #> 5           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #> 6           0|0  0|0           0|0  0|0  0|0           0|0           0|0
    #>   SAMEA12302603 SAMEA12302541 SAMEA12302526 SAMEA12302610 SAMEA12302757
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302542 SAMEA12302711 SAMEA12302720 SAMEA12302672 SAMEA12302753
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302450 SAMEA12302661 SAMEA12302506 SAMEA12302467 SAMEA12302534
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302696 SAMEA12302695 SAMEA12302624 ISIR SAMEA12302670 ISJA
    #> 1           0|0           0|0           0|0  0|0           0|0  0|0
    #> 2           0|0           1|1           0|0  0|0           0|0  0|0
    #> 3           0|0           1|1           0|0  0|0           0|0  0|0
    #> 4           0|0           1|1           0|0  0|0           0|0  0|0
    #> 5           0|0           0|0           0|0  0|0           0|0  0|0
    #> 6           0|0           0|0           0|0  0|0           0|0  0|0
    #>   SAMEA12302679 SAMEA12302538 SAMEA12302507 SAMEA12302587 ISJI SAMEA12302763
    #> 1           0|0           0|0           0|0           0|0  0|0           0|0
    #> 2           0|0           0|0           0|0           0|0  0|0           0|0
    #> 3           0|0           0|0           0|0           0|0  0|0           0|0
    #> 4           0|0           0|0           0|0           0|0  0|0           0|0
    #> 5           0|0           0|0           0|0           0|0  0|0           0|0
    #> 6           0|0           0|0           0|0           0|0  0|0           0|0
    #>   ISJR SAMEA12302478 SAMEA12302787 SAMEA12302818 SAMEA12302837 SAMEA12302799
    #> 1  0|0           0|0           0|0           0|0           0|0           0|0
    #> 2  0|0           0|0           0|0           0|0           0|0           0|0
    #> 3  0|0           0|0           0|0           0|0           0|0           0|0
    #> 4  0|0           0|0           0|0           0|0           0|0           0|0
    #> 5  0|0           0|0           0|0           0|0           0|0           0|0
    #> 6  0|0           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302446 SAMEA12302459 SAMEA12302657 SAMEA12302464 SAMEA12302484
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302584 SAMEA12302570 SAMEA12302578 SAMEA12302573 SAMEA12302468
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302625 SAMEA12302455 SAMEA12302452 SAMEA12302652 SAMEA12302788
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           1|1
    #> 3           0|0           0|0           0|0           0|0           1|1
    #> 4           0|0           0|0           0|0           0|0           1|1
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302487 SAMEA12302514 SAMEA12302653 SAMEA12302691 SAMEA12302676
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           1|1           0|0           0|0
    #> 3           0|0           0|0           1|1           0|0           0|0
    #> 4           0|0           0|0           1|1           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302480 IXYF SAMEA12302752 SAMEA12302516 SAMEA12302828 SAMEA12302518
    #> 1           0|1  0|0           0|0           0|0           0|0           0|0
    #> 2           0|0  0|0           0|0           0|0           0|0           0|0
    #> 3           0|0  0|0           0|0           0|0           0|0           0|0
    #> 4           0|0  0|0           0|0           0|0           0|0           0|0
    #> 5           1|1  0|0           0|0           0|0           0|0           0|0
    #> 6           1|1  0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302814 SAMEA12302616 SAMEA12302586 SAMEA12302795 SAMEA12302826
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           1|1
    #> 3           0|0           0|0           0|0           0|0           1|1
    #> 4           0|0           0|0           0|0           0|0           1|1
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SAMEA12302530 SAMEA12302681 SAMEA12302460 SAMEA12302680 SAMEA12302479
    #> 1           0|0           0|0           0|0           0|0           0|0
    #> 2           0|0           0|0           0|0           0|0           0|0
    #> 3           0|0           0|0           0|0           0|0           0|0
    #> 4           0|0           0|0           0|0           0|0           0|0
    #> 5           0|0           0|0           0|0           0|0           0|0
    #> 6           0|0           0|0           0|0           0|0           0|0
    #>   SRR8474598 IUIW IXPI IXMZ IXYG IXYE SRR8474778 ISKA SRR8474802 IXMJ
    #> 1        0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0
    #> 2        0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0
    #> 3        0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0
    #> 4        0|0  0|0  0|0  0|0  0|0  0|0        0|0  1|1        0|0  0|0
    #> 5        0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0
    #> 6        0|0  0|0  0|0  0|0  0|0  0|0        0|0  0|0        0|0  0|0
    #>   SRR8474801 SRR8474799 SRR8474798 SRR8474797 SRR8474682 SRR8474763 SRR8474679
    #> 1        0|0        0|0        0|0        0|0        0|0        0|0        0|0
    #> 2        1|1        0|0        0|0        0|0        0|0        0|0        0|0
    #> 3        1|1        0|0        0|0        0|0        0|0        0|0        0|0
    #> 4        1|1        0|0        0|0        0|0        0|0        0|0        0|0
    #> 5        0|0        0|0        0|0        0|0        0|0        0|0        0|0
    #> 6        0|0        0|0        0|0        0|0        0|0        0|0        0|0
    #>   SRR8474674 SRR8474675 SRR8474729 SRR8474768 SRR8474761 SRR8474757 SRR8474575
    #> 1        0|0        0|0        0|0        0|0        0|0        0|0        0|0
    #> 2        0|0        0|0        0|0        0|0        1|1        0|0        0|0
    #> 3        0|0        0|0        0|0        0|0        1|1        0|0        0|0
    #> 4        0|0        0|0        0|0        0|0        1|1        0|0        0|0
    #> 5        0|0        0|0        0|0        0|0        0|0        0|0        0|0
    #> 6        0|0        0|0        0|0        0|0        0|0        0|0        0|0
    #>   SRR8474721 IQUS SRR8474691 SRR8474688 SRR8474656 SRR8474652 IQUT JHIM ISGJ
    #> 1        0|0  0|0        0|0        0|0        0|0        0|0  0|0  0|0  0|0
    #> 2        0|0  0|0        0|0        0|0        0|0        0|0  1|1  0|0  0|0
    #> 3        0|0  0|0        0|0        0|0        0|0        0|0  1|1  0|0  0|0
    #> 4        0|0  0|0        0|0        0|0        0|0        0|0  1|1  0|0  0|0
    #> 5        0|0  0|0        0|0        0|0        0|0        0|0  0|0  0|0  0|0
    #> 6        0|0  0|0        0|0        0|0        0|0        0|0  0|0  0|0  0|0
    #>   IQUU IQUW IQUX JHIB JHGW IQUZ IFGS JHFI JHFL JHEU JHHZ ISGS JHFF JHFG JHEQ
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #>   JHFE JHFJ JHFK JHFM SRR8474749 SRR8474747 JHFC SRR8474754 JHER JHET JHFY JHFA
    #> 1  0|0  0|0  0|0  0|0        0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  1|1  0|0        0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  1|1  0|0        0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  1|1  0|0        0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0        0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0        0|0        0|0  0|0        0|0  0|0  0|0  0|0  0|0
    #>   JHFB JHFD SRR8474765 JHGE JHIH JHFU SRR8474640 SRR8474752 JHFW JBPM JHES JHEX
    #> 1  0|0  0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0        0|0  0|0  0|0  0|0        0|0        0|0  0|0  0|0  0|0  0|0
    #>   JHEY JHFX ISIJ JHFT JBPP JHIN JHFH JBPN JHEZ IQWA ISHB IQWB ISMU
    #> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 2  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 3  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 4  0|0  0|0  0|0  0|0  0|0  0|0  1|1  0|0  0|0  0|0  0|0  0|0  0|0
    #> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
    #> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0

``` r
# Example to design a KASP marker on a HIGH impact Deletion variant
library(panGenomeBreedr)
path <- tempdir() # (default directory for saving alignment outputs)

# Path to import sorghum genome sequence for Chromosome 5
path1 <- "https://raw.githubusercontent.com/awkena/panGB/main/Chr05.fa.gz"

# KASP marker design for variant ID: INDEL_Chr05_75106156 in Table 1
lgs1 <- kasp_marker_design(gt_df = gt_region,
                           variant_id_col = 'variant_id',
                           chrom_col = 'chrom',
                           pos_col = 'pos',
                           ref_al_col = 'ref',
                           alt_al_col = 'alt',
                           genome_file = path1,
                           geno_start = 7,
                           marker_ID = "INDEL_Chr05_75106156",
                           chr = "Chr05",
                           plot_draw = TRUE,
                           plot_file = path,
                           region_name = "lgs1")
#> using Gonnet


# View marker alignment output from temp folder
path3 <- file.path(path, list.files(path = path, "alignment_"))
system(paste0('open "', path3, '"')) # Open PDF file from R

on.exit(unlink(path)) # Clear the temp directory on exit
```

The `kasp_marker_design()` function has the following input parameters:

| Argument | Type | Description |
|----|----|----|
| `vcf_file` | `character` | Path to VCF file containing variants. |
| `gt_df` | `data.frame` or `matrix` | Pre-parsed genotype matrix with metadata and sample genotypes. |
| `variant_id_col` | `character` | Column name for variant IDs. |
| `chrom_col` | `character` | Column name for chromosome names. |
| `pos_col` | `character` | Column name for variant positions. |
| `ref_al_col` | `character` | Column name for reference alleles. |
| `alt_al_col` | `character` | Column name for alternate alleles. |
| `geno_start` | `integer` | Column index in `gt_df` where genotype data starts. |
| `marker_ID` | `character` | ID of the variant to design the marker for. |
| `chr` | `character` | Chromosome name used to subset the reference genome. |
| `genome_file` | `character` | Path to genome FASTA file. |
| `plot_draw` | `logical` | Whether to plot upstream/downstream alignment. |
| `plot_file` | `character` | Output path to save the PDF alignment plot. |
| `region_name` | `character` | Optional name for marker region. |
| `maf` | `numeric` | MAF filter for selecting flanking variants. |

The `kasp_marker_design()` function returns a `data.frame` with marker
design metadata:

- `SNP_Name`: Variant ID
- `SNP`: Type of variant (SNP/INDEL)
- `Marker_Name`: Assigned name for the marker
- `Chromosome`: Chromosome name
- `Chromosome_Position`: Variant position
- `Sequence`: Intertek-style polymorphism sequence
- `ReferenceAllele`: Reference allele
- `AlternativeAllele`: Alternate allele

If `plot_draw = TRUE`, a **PDF plot** of sequence alignment will be
saved to `plot_file`.

| <img src='man/figures/alignment.png' align="center" style="width: 700px;" /> |
|:--:|
| *Fig. 2. Alignment of the 100 bp upstream and downstream sequences to the reference genome used for KASP marker design.* |

The required sequence for submission to Intertek for the designed KASP
marker is shown in Table 5.

<table>
<caption>
Table 5: Intertek required sequence for a KASP marker.
</caption>
<thead>
<tr>
<th style="text-align:left;">
SNP_Name
</th>
<th style="text-align:left;">
SNP
</th>
<th style="text-align:left;">
Marker_Name
</th>
<th style="text-align:left;">
Chromosome
</th>
<th style="text-align:right;">
Chromosome_Position
</th>
<th style="text-align:left;">
Sequence
</th>
<th style="text-align:left;">
ReferenceAllele
</th>
<th style="text-align:left;">
AlternativeAllele
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
INDEL_Chr05_75106156
</td>
<td style="text-align:left;">
Deletion
</td>
<td style="text-align:left;">
lgs1
</td>
<td style="text-align:left;">
Chr05
</td>
<td style="text-align:right;">
75106156
</td>
<td style="text-align:left;">
AACCACCGCCACTTCACCGCCGGCGCCGGGGACGACGGCTGCTGGCGCAGAGTAGRAAGTAGRAGTACTCATGCTGCTGCCTGCTGGTGTGTTTCGTATC\[GTAT/-\]GTAGGTGCCTATAAGAGTTTACATGCATGCATCAAAACTAAGGATCCCTCAATGTARYCATGACAGTTCAAAGAAGGCAACAAGGAAGAATTGCTTGCAT
</td>
<td style="text-align:left;">
GTAT
</td>
<td style="text-align:left;">

- </td>
  </tr>
  </tbody>
  </table>

### KASP Marker Validation

The following example demonstrates how to use the customizable functions
in `panGB` to perform hypothesis testing of allelic discrimination for
KASP marker QC and validation.

`panGB` offers customizable functions for KASP marker validation through
hypothesis testing. These functions allow users to easily perform the
following tasks:  
- Import raw or polished KASP genotyping results files (.csv) into R.

- Process imported data and assign FAM and HEX fluorescence colors for
  multiple plates.

- Visualize marker QC using FAM and HEX fluorescence scores for each
  sample.

- Validate the effectiveness of trait-predictive or background markers
  using positive controls.

- Visualize plate design and randomization.

### Reading Raw KASP Full Results Files (.csv)

The `read_kasp_csv()` function allows users to import raw or polished
KASP genotyping full results file (.csv) into R. The function requires
the path of the raw file and the row tags for the different components
of data in the raw file as arguments.

For polished files, the user must extract the `Data` component of the
full results file and save it as a csv file before import.

By default, a typical unedited raw KASP data file uses the following row
tags for genotyping data: `Statistics`, `DNA`, `SNPs`, `Scaling`,
`Data`.

The raw file is imported as a list object in R. Thus, all components in
the imported data can be extracted using the row tag ID as shown in the
code snippet below:

``` r
# Import raw KASP genotyping file (.csv) using the read_kasp_csv() function
library(panGenomeBreedr)

# Set path to the directory where your data is located
# path1 <-  "inst/extdata/Genotyping_141.010_01.csv"
path1 <-  system.file("extdata", "Genotyping_141.010_01.csv",
                       package = "panGenomeBreedr",
                      mustWork = TRUE)

# Import raw data file
file1 <- read_kasp_csv(file = path1, 
                       row_tags = c("Statistics", "DNA", "SNPs", "Scaling", "Data"),
                       data_type = 'raw')

# Get KASP genotyping data for plotting
kasp_dat <- file1$Data
```

### Assigning colors and PCH symbols for KASP cluster plotting

The next step after importing data is to assign FAM and HEX fluorescence
colors to samples based on their observed genotype calls. This step is
accomplished using the `kasp_color()` function in `panGB` as shown in
the code snippet below:

``` r
# Assign KASP fluorescence colors using the kasp_color() function
library(panGenomeBreedr)
# Create a subet variable called plates: masterplate x snpid
  kasp_dat$plates <- paste0(kasp_dat$MasterPlate, '_',
                                 kasp_dat$SNPID)
dat1 <- kasp_color(x = kasp_dat,
                    subset = 'plates',
                    sep = ':',
                    geno_call = 'Call',
                    uncallable = 'Uncallable',
                    unused = '?',
                    blank = 'NTC',
                   assign_cols = c(FAM = "blue", HEX = "gold" , 
                                   het = "forestgreen"))
```

The `kasp_color()` function requires the KASP genotype call file as a
data frame and can do bulk processing if there are multiple master
plates. The default values for the arguments in the `kasp_color()`
function are based on KASP annotations.

The `kasp_color()` function calls the `kasp_pch()` function to
automatically add PCH plotting symbols that can equally be used to group
genotypic clusters on the plot.

When expected genotype calls are available for positive controls in KASP
genotyping samples, we recommend the use of the PCH symbols for grouping
observed genotypes instead of FAM and HEX colors.

The `kasp_color()` function expects that genotype calls are for diploid
state with alleles separated by a symbol. By default KASP data are
separated by `:` symbols.

The `kasp_color()` function returns a list object with the processed
data for each master plate as the components.

### Cluster plot

To test the hypothesis that the designed KASP marker can accurately
discriminate between homozygotes and heterozygotes (allelic
discrimination), a cluster plot needs to be generated.

The `kasp_qc_ggplot()` and `kasp_qc_ggplot2()`functions in `panGB` can
be used to make the cluster plots for each plate and KASP marker as
shown below:

``` r
# KASP QC plot for Plate 05
library(panGenomeBreedr)
kasp_qc_ggplot2(x = dat1[5],
                    pdf = FALSE,
                    Group_id = NULL,
                    scale = TRUE,
                    expand_axis = 0.6,
                    alpha = 0.9,
                    legend.pos.x = 0.6,
                    legend.pos.y = 0.75)
#> $`SE-24-1088_P01_d1_snpSB00804`
```

<div class="figure">

<img src="man/figures/README-plate_05_qc_1-1.png" alt="Fig. 3. Cluster plot for Plate 5 using FAM and HEX colors for grouping observed genotypes." width="100%" />
<p class="caption">
Fig. 3. Cluster plot for Plate 5 using FAM and HEX colors for grouping
observed genotypes.
</p>

</div>

``` r
# KASP QC plot for Plate 05
library(panGenomeBreedr)
 kasp_qc_ggplot2(x = dat1[5],
                  pdf = FALSE,
                  Group_id = 'Group',
                  Group_unknown = '?',
                  scale = TRUE,
                  pred_cols = c('Blank' = 'black', 'False' = 'firebrick3',
                              'True' = 'cornflowerblue', 'Unverified' = 'beige'),
                  expand_axis = 0.6,
                  alpha = 0.9,
                  legend.pos.x = 0.6,
                  legend.pos.y = 0.75)
#> $`SE-24-1088_P01_d1_snpSB00804`
```

<div class="figure">

<img src="man/figures/README-plate_05_qc_2-1.png" alt="Fig. 4. Cluster plot for Plate 5 with an overlay of predictions for positive controls." width="100%" />
<p class="caption">
Fig. 4. Cluster plot for Plate 5 with an overlay of predictions for
positive controls.
</p>

</div>

Color-blind-friendly color combinations are used to visualize verified
genotype predictions (Figure 3).

In Figure 4, the three genotype classes are grouped based on plot PCH
symbols using the FAM and HEX scores for observed genotype calls.

To simplify the verified prediction overlay for the expected genotypes
for positive controls, all possible outcomes are divided into three
categories (TRUE, FALSE, and UNVERIFIED) and color-coded to make it
easier to visualize verified predictions.

BLUE (color code for the TRUE category) means genotype prediction
matches the observed genotype call for the sample.

RED (color code for the FALSE category) means genotype prediction does
not match the observed genotype call for the sample.

BEIGE (color code for the UNVERIFIED category) means three things: an
expected genotype call could not be made before KASP genotyping, or an
observed genotype call could not be made to verify the prediction.

Users can set the `pdf = TRUE` argument to save plots as a PDF file in a
directory outside R. The `kasp_qc_ggplot()` and
`kasp_qc_ggplot2()`functions can generate cluster plots for multiple
plates simultaneously.

To visualize predictions for positive controls to validate KASP markers,
the column name containing expected genotype calls must be provided and
passed to the function using the `Group_id = 'Group'` argument as shown
in the code snippets above. If this information is not available, set
the argument `Group_id = NULL`.

### Summary of Prediction Verification in Plates

The `pred_summary()` function produces a summary of predicted genotypes
for positive controls in each reaction plate after verification (Table
3), as shown in the code snippet below:

``` r
# Get prediction summary for all plates
library(panGenomeBreedr)
my_sum <- pred_summary(x = dat1,
                       snp_id = 'SNPID',
                       Group_id = 'Group',
                       Group_unknown = '?',
                       geno_call = 'Call',
                       rate_out = TRUE)
```

<table>
<caption>
Table 3: Summary of verified prediction status for samples in plates
</caption>
<thead>
<tr>
<th style="text-align:left;">
plate
</th>
<th style="text-align:left;">
snp_id
</th>
<th style="text-align:right;">
false
</th>
<th style="text-align:right;">
true
</th>
<th style="text-align:right;">
unverified
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d1_snpSB00800
</td>
<td style="text-align:left;">
snpSB00800
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.90
</td>
</tr>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d2_snpSB00800
</td>
<td style="text-align:left;">
snpSB00800
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.92
</td>
</tr>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d1_snpSB00803
</td>
<td style="text-align:left;">
snpSB00803
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.66
</td>
</tr>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d2_snpSB00803
</td>
<td style="text-align:left;">
snpSB00803
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.66
</td>
</tr>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d1_snpSB00804
</td>
<td style="text-align:left;">
snpSB00804
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.66
</td>
</tr>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d2_snpSB00804
</td>
<td style="text-align:left;">
snpSB00804
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.66
</td>
</tr>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d1_snpSB00805
</td>
<td style="text-align:left;">
snpSB00805
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.66
</td>
</tr>
<tr>
<td style="text-align:left;">
SE-24-1088_P01_d2_snpSB00805
</td>
<td style="text-align:left;">
snpSB00805
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.66
</td>
</tr>
</tbody>
</table>

The output of the `pred_summary()` function can be visualized as bar
plots using the `pred_summary_plot()` function as shown in the code
snippet below:

``` r
# Get prediction summary for snp:snpSB00804
library(panGenomeBreedr)
my_sum <- my_sum$summ
my_sum <- my_sum[my_sum$snp_id == 'snpSB00804',]

 pred_summary_plot(x = my_sum,
                    pdf = FALSE,
                    pred_cols = c('false' = 'firebrick3', 'true' = 'cornflowerblue',
                                  'unverified' = 'beige'),
                    alpha = 1,
                    text_size = 12,
                    width = 6,
                    height = 6,
                    angle = 45)
#> $snpSB00804
```

<div class="figure">

<img src="man/figures/README-barplot-1.png" alt="Fig. 5. Match/Mismatch rate of predictions for snp: snpSB00804." width="100%" />
<p class="caption">
Fig. 5. Match/Mismatch rate of predictions for snp: snpSB00804.
</p>

</div>

### Plot Plate Design

Users can visualize the observed genotype calls in a plate design format
using the `plot_plate()` function as depicted in Figure 5, using the
code snippet below:

``` r
plot_plate(dat1[5], pdf = FALSE)
#> $`SE-24-1088_P01_d1_snpSB00804`
```

<div class="figure">

<img src="man/figures/README-plate_05_design-1.png" alt="Fig. 6. Observed genotype calls for samples in Plate 5 in a plate design format." width="100%" />
<p class="caption">
Fig. 6. Observed genotype calls for samples in Plate 5 in a plate design
format.
</p>

</div>

# Other Breeder-Centered Functionalities in panGB

`panGB` provides additional functionalities to test hypotheses on the
success of trait introgression pipelines and crosses.

Users can easily generate heatmaps that compare the genetic background
of parents to progenies to ascertain if a target locus was successfully
introgressed or check for the hybridity of F1s. These plots also allow
users to get a visual insight into the amount of parent germplasm
recovered in progenies.

To produce these plots, users must have either polymorphic low or
mid-density marker data and a map file for the markers. **The map file
must contain the marker IDs, their chromosome numbers and positions**.

`panGB`can handle data from KASP, Agriplex and DArTag service providers.

## Working with Agriplex Mid-Density Marker Data

Agriplex data is structurally different from KASP or DArTag data in
terms of genotype call coding and formatting. Agriplex uses `' / '` as a
separator for genotype calls for heterozygotes, and uses single
nucleotides to represent homozygous SNP calls.

## Creating Heatmaps with panGB

To exemplify the steps for creating heatmap, we will use a mid-density
marker data for three groups of near-isogenic lines (NILs) and their
parents (Table 4). The NILs and their parents were genotyped using the
Agriplex platform. Each NIL group was genotyped using 2421 markers.

The imported data frame has the markers as columns and genotyped samples
as rows. It comes with some meta data about the samples. Marker names
are informative: chromosome number and position coordinates are embedded
in the marker names (`Eg. S1_778962: chr = 1, pos = 779862`).

``` r

# Set path to the directory where your data is located
path1 <-  system.file("extdata", "agriplex_dat.csv",
                       package = "panGenomeBreedr",
                      mustWork = TRUE)

# Import raw Agriplex data file
geno <- read.csv(file = path1, header = TRUE, colClasses = c("character")) # genotype calls

library(knitr)
knitr::kable(geno[1:6, 1:10], caption = 'Table 4: Agriplex data format', format = 'html', booktabs = TRUE)
```

<table>
<caption>
Table 4: Agriplex data format
</caption>
<thead>
<tr>
<th style="text-align:left;">
Plate.name
</th>
<th style="text-align:left;">
Well
</th>
<th style="text-align:left;">
Sample_ID
</th>
<th style="text-align:left;">
Batch
</th>
<th style="text-align:left;">
Genotype
</th>
<th style="text-align:left;">
Status
</th>
<th style="text-align:left;">
S1_778962
</th>
<th style="text-align:left;">
S1_1019896
</th>
<th style="text-align:left;">
S1_1613105
</th>
<th style="text-align:left;">
S1_1954298
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RHODES_PLATE1
</td>
<td style="text-align:left;">
D04
</td>
<td style="text-align:left;">
NIL_1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
RTx430a
</td>
<td style="text-align:left;">
Recurrent parent
</td>
<td style="text-align:left;">
A
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
A
</td>
</tr>
<tr>
<td style="text-align:left;">
RHODES_PLATE1
</td>
<td style="text-align:left;">
F04
</td>
<td style="text-align:left;">
NIL_2
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
RTx430b
</td>
<td style="text-align:left;">
Recurrent parent
</td>
<td style="text-align:left;">
A
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
A
</td>
</tr>
<tr>
<td style="text-align:left;">
RHODES_PLATE1
</td>
<td style="text-align:left;">
G04
</td>
<td style="text-align:left;">
NIL_3
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
IRAT204a
</td>
<td style="text-align:left;">
Donor parent
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
A
</td>
</tr>
<tr>
<td style="text-align:left;">
RHODES_PLATE1
</td>
<td style="text-align:left;">
A05
</td>
<td style="text-align:left;">
NIL_4
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
IRAT204b
</td>
<td style="text-align:left;">
Donor Parent
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
A
</td>
</tr>
<tr>
<td style="text-align:left;">
RHODES_PLATE1
</td>
<td style="text-align:left;">
D07
</td>
<td style="text-align:left;">
NIL_5
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
RMES1+\|+\_1
</td>
<td style="text-align:left;">
NIL+
</td>
<td style="text-align:left;">
A
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
A
</td>
</tr>
<tr>
<td style="text-align:left;">
RHODES_PLATE1
</td>
<td style="text-align:left;">
F08
</td>
<td style="text-align:left;">
NIL_6
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
RMES1+\|+\_2
</td>
<td style="text-align:left;">
NIL+
</td>
<td style="text-align:left;">
A
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
A
</td>
</tr>
</tbody>
</table>

To create a heatmap that compares the genetic background of parents and
NILs across all markers, we need to first process the raw Agriplex data
into a numeric format. The panGB package has customizable data wrangling
functions for KASP, Agriplex, and DArTag data.

The `rm_mono()` function can be used to filter out all monomorphic loci
from the data.

Since our imported Agriplex data has informative SNP IDs, we can use the
`parse_marker_ns()` function to generate a map file (Table 5) for the
markers.  
The generated map file is then passed to the `proc_kasp()` function to
order the SNP markers according to their chromosome numbers and
positions.

The `kasp_numeric()` function converts the output of the `proc_kasp()`
function into a numeric format (Table 6). The re-coding to numeric
format is done as follows:

- Homozygous for Parent 1 allele = 1.
- Homozygous for Parent 2 allele = 0.
- Heterozygous = 0.5.
- Monomorphic loci = -1.
- Loci with a suspected genotype error = -2.
- Loci with at least one missing parental or any other genotype = -5.

``` r

# Parse snp ids to generate a map file
library(panGenomeBreedr)

# Data for stg5 NILs
stg5 <- geno[geno$Batch == 3, -c(1:6)] 
rownames(stg5) <- geno$Genotype[17:25]

# Remove monomorphic loci from data
stg5 <- rm_mono(stg5)

# Parse snp ids to generate a map file
snps <- colnames(stg5) # Get snp ids
map_file <- parse_marker_ns(x = snps, sep = '_', prefix = 'S')

# order markers in map file
map_file <- order_markers(x = map_file)
```

<table>
<caption>
Table 5: Map file for the imported Agriplex data.
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
snpid
</th>
<th style="text-align:right;">
chr
</th>
<th style="text-align:right;">
pos
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1.317
</td>
<td style="text-align:left;">
S1_402592
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
402592
</td>
</tr>
<tr>
<td style="text-align:left;">
1.1
</td>
<td style="text-align:left;">
S1_778962
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
778962
</td>
</tr>
<tr>
<td style="text-align:left;">
1.633
</td>
<td style="text-align:left;">
S1_825853
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
825853
</td>
</tr>
<tr>
<td style="text-align:left;">
1.318
</td>
<td style="text-align:left;">
S1_1218846
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1218846
</td>
</tr>
<tr>
<td style="text-align:left;">
1.2
</td>
<td style="text-align:left;">
S1_1613105
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1613105
</td>
</tr>
<tr>
<td style="text-align:left;">
1.319
</td>
<td style="text-align:left;">
S1_1727150
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1727150
</td>
</tr>
<tr>
<td style="text-align:left;">
1.3
</td>
<td style="text-align:left;">
S1_1954298
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1954298
</td>
</tr>
<tr>
<td style="text-align:left;">
1.4
</td>
<td style="text-align:left;">
S1_1985365
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1985365
</td>
</tr>
</tbody>
</table>

``` r
# Process genotype data to re-order SNPs based on chromosome and positions
stg5 <- proc_kasp(x = stg5,
                  kasp_map = map_file,
                  map_snp_id = "snpid",
                  sample_id = "Genotype",
                  marker_start = 1,
                  chr = 'chr',
                  chr_pos = 'pos')

# Convert to numeric format for plotting
num_geno <- kasp_numeric(x = stg5,
                         rp_row = 1,
                         dp_row = 3,
                         sep = ' / ',
                         data_type = 'agriplex')
```

<table>
<caption>
Table 6: Agriplex data converted to a numeric format.
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
S1_402592
</th>
<th style="text-align:right;">
S1_778962
</th>
<th style="text-align:right;">
S1_825853
</th>
<th style="text-align:right;">
S1_1218846
</th>
<th style="text-align:right;">
S1_1613105
</th>
<th style="text-align:right;">
S1_1727150
</th>
<th style="text-align:right;">
S1_1954298
</th>
<th style="text-align:right;">
S1_1985365
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BTx623a
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
BTx623b
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
BTx642a
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
BTx642b
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Stg5+\|+\_1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Stg5+\|+\_2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Stg5-\|-\_1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Stg5-\|-\_2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Stg5-\|-\_3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

All is now set to generate the heatmap (Figure 6) using the
`cross_qc_ggplot()` function, as shown in the code snippet below:

``` r

# Get prediction summary for snp:snpSB00804
library(panGenomeBreedr)
# Create a heatmap that compares the parents to progenies
cross_qc_ggplot(x = num_geno,
                map_file = map_file,
                snp_ids = 'snpid',
                chr = 'chr',
                chr_pos = 'pos',
                parents = c("BTx623a", "BTx642a"),
                pdf = FALSE,
                filename = 'background_heatmap',
                legend_title = 'stg5_NILs',
                alpha = 0.9,
                text_size = 15)
#> $Batch1
```

<div class="figure">

<img src="man/figures/README-heatmap1-1.png" alt="Fig. 6. A heatmap that compares the genetic background of parents and stg5 NIL progenies across all markers." width="100%" />
<p class="caption">
Fig. 6. A heatmap that compares the genetic background of parents and
stg5 NIL progenies across all markers.
</p>

</div>

The `cross_qc_ggplot()` function is a wrapper for functions in the
`ggplot2` package.

Users must specify the IDs for the two parents using the `parents`
argument. In the code snippet above, the recurrent parent is `BTx623`
and the donor parent for the *stg5* locus is `BTx642`.

The `group_sz` argument must be specified to plot the heatmap in batches
of progenies to avoid cluttering the plot with many observations.

Users can set the `pdf = TRUE` argument to save plots as a PDF file in a
directory outside R.

## Trait Introgression Hypothesis Testing

To test the hypothesis that the *stg5* NIL development was effective, we
can use the `cross_qc_annotate()` function to generate a heatmap (Figure
7) with an annotation of the position of the *stg5* locus on Chr 1, as
shown below:

``` r

###########################################################################
# Subset data for the first 30 markers on Chr 1
stg5_ch1 <- num_geno[, map_file$chr == 1][,1:30] 

# Get the map file for subset data
stg5_ch1_map <- parse_marker_ns(colnames(stg5_ch1))

# Annotate a heatmap to show the stg5 locus on Chr 1
# The locus is between positions 0.98 - 1.8 Mbp on Chr 1
cross_qc_annotate(x = stg5_ch1,
                  map_file = stg5_ch1_map,
                  snp_ids = 'snpid',
                  chr = 'chr',
                  chr_pos = 'pos',
                  parents = c("BTx623a", "BTx642a"),
                  trait_pos = list(stg5 = c(start = .98e6, end = 1.8e6)),
                  text_scale_fct = 0.3,
                  pdf = FALSE,
                  legend_title = 'Stg5_NILs',
                  alpha = 0.9,
                  text_size = 15)
#> $Batch1
```

<div class="figure">

<img src="man/figures/README-heatmap2-1.png" alt="Fig. 7. Heatmap annotation of the stg5 locus on Chr 1." width="100%" />
<p class="caption">
Fig. 7. Heatmap annotation of the stg5 locus on Chr 1.
</p>

</div>

In the code snippet above, the numeric matrix of genotype calls and its
associated map file are required.

The recurrent and donor parents must be specified using the `parents`
argument.

The `snp_ids, chr, and chr_pos` arguments can be used to specify the
column names for marker IDs, chromosome number and positions in the
attached map file.  
The `trait_pos` argument was used to specify the position of the target
locus (*stg5*) on chromosome one. Users can specify the positions of
multiple target loci as components of a list object for annotation.

In Figure 7, the color intensity correlates positively with the marker
density or coverage. Thus, areas with no color (white vertical gaps)
depicts gaps in the marker coverage in the data.

## Decision Support for Marker-Assisted Backcrossing in panGB

Users can use the `calc_rpp_bc()` function in `panGB` to calculate the
proportion of recurrent parent background (RPP) fully recovered in
backcross progenies.

In the computation, partially regions are ignored, hence, heterozygous
scores are not used.

The output for he `calc_rpp_bc()` function can be passed to the
`rpp_barplot()` function to visualize the computed RPP values for
progenies as a bar plot. Users can specify an RPP threshold to easily
identify lines that have RPP values above or equal to the defined RPP
threshold on the bar plot.

We can compute and visualize the observed RPP values for the *stg5* NILs
across all polymorphic loci as shown in the code snippet below:

``` r

# Calculate weighted RPP
rpp <- calc_rpp_bc(x = num_geno,
                   map_file = map_file,
                   map_chr = 'chr',
                   map_pos = 'pos',
                   map_snp_ids = 'snpid',
                   rp = 1,
                   rp_num_code = 1,
                   na_code = -5,
                   weighted = TRUE)

# Generate bar plot for RPP values
rpp_barplot(rpp_df = rpp,
            rpp_threshold = 0.93,
            text_size = 18,
            text_scale_fct = 0.1,
            alpha = 0.9,
            bar_width = 0.5,
            aspect_ratio = 0.5,
            pdf = FALSE)
```

<div class="figure">

<img src="man/figures/README-barplot_rpp1-1.png" alt="Fig. 8. Computed RPP values for the stg5 NILs." width="100%" />
<p class="caption">
Fig. 8. Computed RPP values for the stg5 NILs.
</p>

</div>

The `calc_rpp_bc()` function in `panGB` provides two algorithms for
computing the observed RPP values: weighted and unweighted RPP values.
We recommend the use of the weighted algorithm to account for
differences in the marker coverage across the genome.

The algorithm for the weighted RPP values is explained below.

### Weighted RPP Computation in panGB

Let $w_i$ represent the weight for marker $i$, based on the relative
distances to its adjacent markers.

For a set of markers with positions $p_1, p_2, \ldots, p_n$, where
$d_i = p_{i+1} - p_i$ represents the distance between adjacent markers,
the weights can be calculated as follows:

1.  **For the first marker** $i = 1$:

    $$w_1 = \frac{d_1}{2 \sum_{i=1}^{n-1} d_i}$$

2.  **For a middle marker** $1 < i < n$:

    $$w_i = \frac{d_{i-1} + d_i}{2 \sum_{i=1}^{n-1} d_i}$$

3.  **For the last marker** $i = n$:

    $$w_n = \frac{d_{n-1}}{2 \sum_{i=1}^{n-1} d_i}$$

where:

- $d_i$ is the distance between marker $i$ and marker $i+1$,
- $sum_{i=1}^{n-1} d_i$ is the total distance across all segments, used
  for normalization.

Let $RPP$ represent the Recurrent Parent Proportion based on relative
distance weighting. If $w_i$ is the weight for each marker $i$, and
$m_i$ represents whether marker $i$ matches the recurrent parent
$m_i = 1$ if it matches, $m_i = 0$ otherwise), then the weighted RPP is
calculated as:

$$RPP_{weighted} = \sum_{i=1}^n w_i\cdot m_i$$

The unweighted RPP is calculated without the use of the weights as
follows:

$$RPP_{unweighted} = \frac{\sum_{i=1}^n m_i} n$$

where:

- $w_i$ is the weight of marker $i$, calculated based on the relative
  distance it covers,
- $m_i$ is the match indicator for marker $i$ (1 if matching the
  recurrent parent, 0 otherwise),
- $n$ is the total number of markers.

This formula provides the sum of the weighted contributions from each
marker, representing the proportion of the recurrent parent genome in
the individual.

## Troubleshooting

If the package does not run as expected, check the following:

- Was the package properly installed?

- Do you have the required dependencies installed?

- Were any warnings or error messages returned during package
  installation?

- Are all packages up to date before installing panGB?

# Authors and contributors

- [Alexander Wireko Kena](https://www.github.com/awkena)

- [Cruet Burgos](https://www.morrislab.org/people/clara-cruet-burgos)

- [Linly Banda](https://www.biofortificationlab.org/people/linly-banda)

- [Jacques
  Faye](https://sites.google.com/site/morrislaboratory/people/jacques-faye)

- [Fanna Maina](https://www.morrislab.org/people/fanna-maina)

- [Terry
  Felderhoff](https://www.agronomy.k-state.edu/about/people/faculty/felderhoff-terry/)

- [Geoffrey Preston
  Morris](https://www.morrislab.org/people/geoff-morris)

# License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

# Support and Feedback

For support and submission of feedback, email the maintainer **Alexander
Kena, PhD** at <alex.kena24@gmail.com>
