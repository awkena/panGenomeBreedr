# List all tables in the database

This function lists all available tables from either a local database
connection or the online database, depending on the selected mode.

## Usage

``` r
list_tables(con = NULL, connect_db_mode = c("local", "online"))
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A character vector of table names in the database.

## Examples

``` r
# \donttest{
# --- Online Mode ---
# Load the package
library(panGenomeBreedr)

# List all tables in the online database
online_tables <- list_tables(connect_db_mode = 'online')
print(online_tables)
#> [1] "annotations" "genotypes"   "metadata"    "variants"   

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# List tables from the local connection (default mode is 'local')
local_tables <- list_tables(con = con_local)
# print(local_tables)

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
