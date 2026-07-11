# Compute allele frequencies for a genotype matrix (online)

Compute allele frequencies for a genotype matrix (online)

## Usage

``` r
pg_calc_af(gt, variant_id_col = "variant_id", chrom_col = NULL, pos_col = NULL)
```

## Arguments

- gt:

  A data frame or matrix where rows are variants and columns are
  samples. May also contain metadata columns (ID, Chrom, Pos).

- variant_id_col:

  Character. Name of the column containing variant IDs. Default is
  'variant_id'.

- chrom_col:

  Character. Optional name of the chromosome column.

- pos_col:

  Character. Optional name of the position column.

## Value

A data frame containing the variant metadata and two calculated columns:

- `ref_af`: Reference allele frequency.

- `alt_af`: Alternate allele frequency.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# First, query genotype data for a region from the online database
# This will be the 'gt' (genotype matrix) input for pg_calc_af
gt_data_region <- pg_query_db(
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)

# Calculate allele frequencies for the fetched genotype data
af_metrics <- pg_calc_af(
  gt = gt_data_region,
  variant_id_col = "variant_id",
  chrom_col = "chrom",
  pos_col = "pos"
)
print(head(af_metrics))
} # }
```
