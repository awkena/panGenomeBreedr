# Extract variants based on functional mutation impact (online)

Extract variants based on functional mutation impact (online)

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

  A character vector specifying the variant impact types to retain.
  Allowed classifications are `"HIGH"`, `"MODERATE"`, `"LOW"`, and
  `"MODIFIER"`.

- chrom:

  A character value specifying the chromosome name (e.g., `"Chr05"`).

- start:

  Integer. Start coordinate for the genomic region.

- end:

  Integer. End coordinate for the genomic region.

## Value

A data frame containing variant information and associated functional
annotations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Query for high-impact variants in a specific genomic region
high_impact_vars <- pg_query_by_impact(
  impact_level = "HIGH",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)
print(high_impact_vars)

} # }
```
