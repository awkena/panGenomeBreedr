# Extract variants based on allele frequencies within a genomic region

This function connects to the public panGenomeBreedr API to query
genotypes within a specific genomic range and filters the results to
only include variants within the specified alternate allele frequency
(AF) thresholds.

## Usage

``` r
pg_query_by_af(min_af = 0, max_af = 1, chrom = NULL, start = NULL, end = NULL)
```

## Arguments

- min_af:

  Numeric. Minimum alternate allele frequency (0-1). Default is 0.

- max_af:

  Numeric. Maximum alternate allele frequency (0-1). Default is 1.

- chrom:

  Character. Chromosome name (e.g., "Chr05").

- start, end:

  Numeric. Genomic start and end coordinates.

## Value

A data frame containing variant metadata (ID, Chrom, Pos) and the
calculated `ref_af` and `alt_af`.
