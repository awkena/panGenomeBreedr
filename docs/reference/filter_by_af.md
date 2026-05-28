# Filter extracted genotypes based on alternate allele frequency (local)

Passes an extracted genotype matrix through calculation utilities to
compute population-level allele frequencies, then subsets rows down to
variants strictly falling between the specified minimum and maximum
alternate allele frequency (MAF) thresholds.

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

  A matrix or data frame containing variants (rows) across samples
  (columns), including required variant identification and optional
  genomic coordinate metadata columns.

- variant_id_col:

  Character. The column name in `gt` matching unique variant
  identifiers. Defaults to `"variant_id"`.

- chrom_col:

  Character. The column name in `gt` matching chromosome names. Defaults
  to `"chrom"`.

- pos_col:

  Character. The column name in `gt` matching genomic positions.
  Defaults to `"pos"`.

- min_af:

  Numeric. Minimum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `0`.

- max_af:

  Numeric. Maximum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `1`.

## Value

A data frame of variants passing the frequency thresholds, tracking
their coordinates, calculated reference allele frequencies (`ref_af`),
and alternate allele frequencies (`alt_af`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Connect to the database pipeline
con <- connect_local_db(folder_path = "~/Desktop/curated_sorghum_variant_resource")

# Pull wide genotype records for a locus window range
gt_region <- query_db(con = con, chrom = "Chr05", start = 75104537, 
                      end = 75106403, table_name = "genotypes")

# Filter out rare variants or monomorphic lines by restricting the alternate allele frequency window
filtered_vars <- filter_by_af(gt_region, variant_id_col = "variant_id", 
                              chrom_col = "chrom", pos_col = "pos", 
                              min_af = 0.05, max_af = 0.95)

# Safely close out the database connection
disconnect_local_db(con)
} # }
```
