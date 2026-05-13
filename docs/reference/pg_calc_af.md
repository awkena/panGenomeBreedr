# Compute allele frequencies for a genotype matrix

This function calculates the reference and alternate allele frequencies
for a variants-by-samples matrix. It handles both phased (\|) and
unphased (/) VCF genotype strings.

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
