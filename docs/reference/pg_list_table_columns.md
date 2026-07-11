# Get table schema details (online)

Get table schema details (online)

## Usage

``` r
pg_list_table_columns(
  table_name = c("variants", "annotations", "genotypes", "metadata")
)
```

## Arguments

- table_name:

  A character value specifying the name of the table. Defaults to
  "variants". Available options are "variants", "annotations",
  "genotypes", or "metadata".

## Value

A data frame containing column metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Get column metadata for the "variants" table
pg_list_table_columns_result <- pg_list_table_columns(table_name = "variants")
print(pg_list_table_columns_result)
} # }
```
