# Get sample metadata from the database (online)

Get sample metadata from the database (online)

## Usage

``` r
pg_get_sample_metadata(query_col = NULL, query_value = NULL)
```

## Arguments

- query_col:

  Character. Optional metadata column name to filter by (e.g.,
  `"countryorigin"`). If `NULL`, all database metadata records are
  returned.

- query_value:

  Character or Numeric. The explicit attribute value to match within
  `query_col`.

## Value

A data frame containing the sample metadata records.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Retrieve all sample metadata from the online database
all_metadata <- pg_get_sample_metadata()
print(head(all_metadata))

# Retrieve sample metadata filtered by 'countryorigin' = "Ghana"
ghana_metadata <- pg_get_sample_metadata(
  query_col = "countryorigin",
  query_value = "Ghana"
)
print(head(ghana_metadata))
} # }

```
