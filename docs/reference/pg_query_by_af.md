# Query variants by allele frequency (online)

Query variants by allele frequency (online)

## Usage

``` r
pg_query_by_af(min_af = 0, max_af = 1, chrom = NULL, start = NULL, end = NULL)
```

## Arguments

- min_af:

  Numeric. Minimum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `0`.

- max_af:

  Numeric. Maximum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `1`.

- chrom:

  A character value specifying the target chromosome name (e.g.,
  `"Chr05"`).

- start:

  Integer. Start coordinate for the target window region.

- end:

  Integer. End coordinate for the target window region.

## Value

A data frame containing variant metadata and the calculated allele
frequencies, filtered by the specified thresholds.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Query variants from the online database for a specific region
# and filter for common variants (MAF >= 5%)
common_variants_online <- pg_query_by_af(
  min_af = 0.05,
  max_af = 0.95,
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)

# Print the head of the resulting data frame
print(head(common_variants_online))
} # }
```
