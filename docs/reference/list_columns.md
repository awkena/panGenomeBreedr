# Get table schema details

Retrieves the schema (column names and data types) for a specified
table. It can operate in 'local' or 'online' mode.

## Usage

``` r
list_columns(
  con = NULL,
  table_name = c("variants", "annotations", "genotypes", "metadata"),
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- table_name:

  A character value specifying the name of the table to inspect. Must be
  one of \`"variants"\`, \`"annotations"\`, \`"genotypes"\`, or
  \`"metadata"\`.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame containing the schema details of the target table.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.

# Inspect the schema of the local "variants" table
local_schema <- list_columns(con = con_local, table_name = "variants")
print(local_schema)
#>    column_name column_type null  key default extra
#> 1   variant_id     VARCHAR  YES <NA>    <NA>  <NA>
#> 2        chrom     VARCHAR  YES <NA>    <NA>  <NA>
#> 3          pos      DOUBLE  YES <NA>    <NA>  <NA>
#> 4          ref     VARCHAR  YES <NA>    <NA>  <NA>
#> 5          alt     VARCHAR  YES <NA>    <NA>  <NA>
#> 6         qual     VARCHAR  YES <NA>    <NA>  <NA>
#> 7       filter     VARCHAR  YES <NA>    <NA>  <NA>
#> 8 variant_type     VARCHAR  YES <NA>    <NA>  <NA>

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To query the public online resource instead:
online_schema <- list_columns(table_name = "variants", connect_db_mode = 'online')
print(online_schema)
} # }
```
