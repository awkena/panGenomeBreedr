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
# Load the package
library(panGenomeBreedr)

# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.

# Get row counts for all mounted Parquet views
local_summary <- summarize_database(con = con_local)
print(local_summary)
#>                   table n_rows
#> 1           annotations    421
#> 2             genotypes    130
#> 3              metadata   1676
#> 4            pcil_genes      2
#> 5       pcil_genomewide     47
#> 6              pcil_ibs   1081
#> 7   pcil_introgressions    473
#> 8         pcil_metadata   1383
#> 9  pcil_sample_metadata     47
#> 10             variants    130

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To query the public online resource instead:
online_summary <- summarize_database(connect_db_mode = 'online')
print(online_summary)
} # }
```
