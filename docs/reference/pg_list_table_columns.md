# Check the column names and types for any table in the pangenome database.

This function connects to the public panGenomeBreedr API and retrieves
metadata about the columns in a specified table.

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
