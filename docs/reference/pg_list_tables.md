# List all tables in the database (online)

List all tables in the database (online)

## Usage

``` r
pg_list_tables()
```

## Value

A character vector of table names.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# List all tables in the online database
pg_list_tables_result <- pg_list_tables()
print(pg_list_tables_result)
#> [1] "genotypes"   "variants"    "annotations" "metadata"   
# }
```
