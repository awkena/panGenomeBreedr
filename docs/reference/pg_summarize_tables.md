# Get table names and row counts from the database (online)

Get table names and row counts from the database (online)

## Usage

``` r
pg_summarize_tables()
```

## Value

A data frame with two columns: `table` and `n_rows`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Summarize all tables in the online database
pg_summarize_tables_result <- pg_summarize_tables()
print(pg_summarize_tables_result)
} # }
```
