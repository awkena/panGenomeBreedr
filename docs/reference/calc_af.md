# Compute allele frequencies for a VCF genotype matrix (variant x samples). Chromosome and position may be included in the data.

Compute allele frequencies for a VCF genotype matrix (variant x
samples). Chromosome and position may be included in the data.

## Usage

``` r
calc_af(gt, variant_id_col = "variant_id", chrom_col = NULL, pos_col = NULL)
```

## Arguments

- gt:

  A matrix or data frame of variants x samples with chromosome and
  position meta data for each variant.

- variant_id_col, chrom_col, pos_col:

  A character value specifying the column names of variant IDs,
  chromosome, and positions in `gt`.

## Value

A data frame of variants with their reference and alternate allele
frequencies.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)

# Define tempdir
path <- tempdir()

# Mini SQLite database
mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
                     package = "panGenomeBreedr",
                     mustWork = TRUE)

# Path to SQLite databases: INDEL and SNP
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')

# Unzip compressed mini database and save in tempdir
R.utils::gunzip(mini_db,
               destname = mini_db_path,
               remove = FALSE)

# Extract variants within a candidate gene and their calculate allele
# frequencies
variant_region <- query_db(db_path = mini_db_path,
                              chrom = "Chr05",
                              start = 75104537,
                              end = 75106403,
                              table_name = "genotypes",
                              gene_name = "Sobic.005G213600") |>
calc_af(variant_id_col = 'variant_id',
       chrom_col = 'chrom',
       pos_col = 'pos')

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
