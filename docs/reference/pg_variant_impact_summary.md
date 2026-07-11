# Get summary of database variant impacts by chromosome (online)

Get summary of database variant impacts by chromosome (online)

## Usage

``` r
pg_variant_impact_summary()
```

## Value

A data frame in wide format where each row is a chromosome and columns
represent the counts for each impact category.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Get variant impact summary from the online database
pg_variant_impact_summary_result <- pg_variant_impact_summary()
print(pg_variant_impact_summary_result)
} # }
```
