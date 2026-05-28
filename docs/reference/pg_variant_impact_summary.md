# Summarize variant mutation impacts per chromosome (online)

This function connects to the public panGenomeBreedr API and summarizes
the distribution of mutation impacts across chromosomes.

## Usage

``` r
pg_variant_impact_summary()
```

## Value

A data frame in wide format where each row is a chromosome and columns
represent the counts for each impact category.
