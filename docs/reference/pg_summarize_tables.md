# Name and row count for each table in the pangenome database.

This function connects to the public panGenomeBreedr API and returns a
summary data frame containing the table names and their respective row
counts.

## Usage

``` r
pg_summarize_tables()
```

## Value

A data frame with two columns: `table` and `n_rows`.
