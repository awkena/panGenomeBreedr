# Filter extracted variants based on alternate allele frequency.

Filter extracted variants based on alternate allele frequency.

## Usage

``` r
filter_by_af(
  gt,
  variant_id_col = "variant_id",
  chrom_col = "chrom",
  pos_col = "pos",
  min_af = 0,
  max_af = 1
)
```

## Arguments

- gt:

  A matrix or data frame of variants x samples with chromosome and
  position meta data for each variant.

- variant_id_col, chrom_col, pos_col:

  A character value specifying the column names of variant IDs,
  chromosome, and positions in `gt`.

- min_af, max_af:

  A numeric value indicating the minimum and maximum allele frequency
  thresholds for filtering variants.

## Value

A data frame of filtered variants with their reference and alternate
allele frequencies.

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

# Extract variants within a candidate gene and filter by allele frequency
variant_region <- query_db(db_path = mini_db_path,
                              chrom = "Chr05",
                              start = 75104537,
                              end = 75106403,
                              table_name = "genotypes",
                              gene_name = "Sobic.005G213600") |>
filter_by_af(variant_id_col = 'variant_id',
             chrom_col = 'chrom',
             pos_col = 'pos',
             min_af = 0.01,
             max_af = 0.99)

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
