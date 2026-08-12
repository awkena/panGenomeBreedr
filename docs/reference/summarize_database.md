# Get table names and row counts from the database

Queries the active database to retrieve a list of all mounted tables
along with their respective total row counts. It can operate in 'local'
or 'online' mode.

## Usage

``` r
summarize_database(con = NULL, connect_db_mode = c("local", "online"))
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame containing two columns: `table` (the table name) and
`n_rows` (the total number of rows in that table).

## Examples

``` r
# \donttest{
# --- Online Mode ---
# Load the package
library(panGenomeBreedr)

# Summarize all tables in the online database
online_summary <- summarize_database(connect_db_mode = 'online')
print(online_summary)
#>         table   n_rows
#> 1 annotations 61679112
#> 2    metadata     1676
#> 3    variants 38230911
#> 4   genotypes 38230930

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Get row counts for all mounted Parquet views
local_summary <- summarize_database(con = con_local)
# print(local_summary)

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
