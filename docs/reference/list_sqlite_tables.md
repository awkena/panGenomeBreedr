# List all tables in the database (local)

List all tables in the database (local)

## Usage

``` r
list_sqlite_tables(con)
```

## Arguments

- con:

  A DBI connection object to the local database.

## Value

A character vector of table names in the database.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

my_db_folder <- "~/Desktop/curated_sorghum_variant_resource"
con <- connect_local_db(folder_path = my_db_folder)
list_sqlite_tables(con)
disconnect_local_db(con)
} # }
```
