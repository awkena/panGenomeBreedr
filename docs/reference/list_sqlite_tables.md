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

# Define the path to the curated_sorghum_variant_resource folder
my_db_folder <- "~/Desktop/curated_sorghum_variant_resource"

# Establish the connection
con <- connect_local_db(folder_path = my_db_folder)

# List tables
list_sqlite_tables_result <- list_sqlite_tables(con)
print(list_sqlite_tables_result)

# Disconnect at the end of your session
disconnect_local_db(con)

} # }
```
