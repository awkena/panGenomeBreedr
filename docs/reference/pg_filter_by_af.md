# Filter extracted variants based on alternate allele frequency

Calculates allele frequencies for a genotype matrix and filters variants
based on a user-defined range locally on your machine. Useful for
removing monomorphic or rare variants from pangenome queries.

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

  A data frame or matrix of variants x samples, typically the output
  from
  [`pg_query_db`](https://awkena.github.io/panGenomeBreedr/reference/pg_query_db.md).

- variant_id_col:

  Character. Column name for variant IDs. Default is 'variant_id'.

- chrom_col:

  Character. Optional column name for chromosome.

- pos_col:

  Character. Optional column name for genomic position.

- min_af:

  Numeric. Minimum alternate allele frequency threshold (0-1).

- max_af:

  Numeric. Maximum alternate allele frequency threshold (0-1).

## Value

A data frame containing variant metadata and calculated frequencies for
variants that passed the filter.

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Query region via the API and pipe into filter to remove rare variants (MAF < 0.05)
# Notice: No database connection needed!
filtered_vars <- pg_query_db(
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
) |>
  pg_filter_by_af(min_af = 0.05, max_af = 0.95)
} # }
```
