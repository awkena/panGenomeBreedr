# Filter genotypes by allele frequency (online)

Filter genotypes by allele frequency (online)

## Usage

``` r
pg_filter_by_af(
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

  A data frame containing genotype calls queried from the genotypes
  table.

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

A data frame containing variant metadata and calculated frequencies for
variants that passed the filter.

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Query genotypes table and filter to remove rare variants (MAF < 0.05)
filtered_vars <- pg_query_db(
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
) |>
  pg_filter_by_af(min_af = 0.05, max_af = 0.95)
} # }
```
