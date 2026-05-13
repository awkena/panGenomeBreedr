# Get variant statistics from the PostgreSQL database

This function connects to the public panGenomeBreedr API and calculates
summary statistics for variants per chromosome, including variant counts
and genomic ranges.

## Usage

``` r
pg_variant_stats(include_annotations = TRUE)
```

## Arguments

- include_annotations:

  A logical value indicating whether to include statistics for the
  annotations table. Defaults to `TRUE`.

## Value

A data frame containing variant statistics.
