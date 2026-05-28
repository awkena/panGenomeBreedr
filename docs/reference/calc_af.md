# Compute allele frequencies for a genotype matrix (local)

Calculates the reference (REF) and alternate (ALT) allele frequencies
for an extracted wide matrix or data frame of variants by samples.
Converts standard VCF genotype strings (e.g., "0/0", "0\|1", "1\|1")
into numerical alternate allele dosages (0, 1, 2) to perform
population-level frequency calculations.

## Usage

``` r
calc_af(gt, variant_id_col = "variant_id", chrom_col = NULL, pos_col = NULL)
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

  Character. Optional column name in `gt` matching chromosome names.
  Defaults to `NULL`.

- pos_col:

  Character. Optional column name in `gt` matching genomic positions.
  Defaults to `NULL`.

## Value

A data frame containing unique variant identifiers, optional coordinate
tracking columns, calculated reference allele frequencies (`ref_af`),
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

# Compute allele frequencies across all samples in the dataset
af_metrics <- calc_af(gt_region, variant_id_col = "variant_id", 
                      chrom_col = "chrom", pos_col = "pos")

# Safely close out the database connection
disconnect_local_db(con)
} # }
```
