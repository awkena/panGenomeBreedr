# Extract variants from the annotation table based on impact type.

This function connects to the public panGenomeBreedr API to retrieve
variants filtered by snpEff impact levels (HIGH, MODERATE, etc.) and
optional genomic coordinates.

## Usage

``` r
pg_query_by_impact(
  impact_level = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
  chrom = NULL,
  start = NULL,
  end = NULL
)
```

## Arguments

- impact_level:

  Character vector. One or more of "HIGH", "MODERATE", "LOW",
  "MODIFIER". Defaults to all four.

- chrom:

  Character. Optional chromosome name (e.g., "Chr05").

- start, end:

  Numeric. Optional genomic start and end positions.

## Value

A data frame containing variant information and associated functional
annotations.
