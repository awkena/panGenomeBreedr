# Get variant statistics from the database (online)

Get variant statistics from the database (online)

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Get variant statistics from the online database, including annotation counts
pg_variant_stats_result <- pg_variant_stats(include_annotations = TRUE)
print(pg_variant_stats_result)

# Get variant statistics without annotation counts
pg_variant_stats_no_ann <- pg_variant_stats(include_annotations = FALSE)
print(pg_variant_stats_no_ann)
} # }
```
