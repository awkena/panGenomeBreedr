# Compute allele frequencies for a genotype matrix (local)

Compute allele frequencies for a genotype matrix (local)

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

# Define the path to the curated_sorghum_variant_resource folder
my_db_folder <- "~/Desktop/curated_sorghum_variant_resource"

# Establish the connection
con <- connect_local_db(folder_path = my_db_folder)

# Pull wide genotype records for a locus window range
gt_region <- query_db(con = con, chrom = "Chr05", start = 75104537, 
                      end = 75106403, table_name = "genotypes")

# Compute allele frequencies across all samples in the dataset
af_metrics <- calc_af(gt_region, variant_id_col = "variant_id", 
                      chrom_col = "chrom", pos_col = "pos")

# Disconnect at the end of your session
disconnect_local_db(con)
} # }
```
