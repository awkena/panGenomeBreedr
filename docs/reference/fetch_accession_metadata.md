# Get sample metadata from the database

Retrieves sample metadata for the genotypes table from the database. It
can operate in 'local' mode to query a local database connection, or
'online' mode to fetch data from the API. The function can return all
metadata or filter by a specific column-value pair.

## Usage

``` r
fetch_accession_metadata(
  con = NULL,
  query_col = NULL,
  query_value = NULL,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- query_col:

  Character. Optional metadata column name to filter by (e.g.,
  \`"countryorigin"\`). If \`NULL\`, all metadata records are returned.

- query_value:

  Character or Numeric. The value to match within \`query_col\`.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame containing the sample metadata records.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# --- Online Mode ---
# Retrieve all sample metadata from the online database
all_metadata_online <- fetch_accession_metadata(connect_db_mode = 'online')
print(all_metadata_online[1:3, 1:5])
#>   array_index  lib         sample pinumber isnum
#> 1           1 IDMM           Juar     <NA>  <NA>
#> 2           2 ISGC     00KOF5DT19     <NA>  <NA>
#> 3           3 ISGK 02-SB-F4DT-275     <NA>  <NA>

# Retrieve metadata for samples from Ghana from the online database
ghana_metadata_online <- fetch_accession_metadata(
  query_col = "countryorigin",
  query_value = "Ghana",
  connect_db_mode = 'online'
)
print(ghana_metadata_online[1:3, 1:5])
#>   array_index  lib   sample pinumber   isnum
#> 1         486 IQST IS_25055 Grif9991 IS25055
#> 2        1234 ISMA PI585448 PI585448 IS25050
#> 3        1235 ISMG PI585452 PI585452 IS25059

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db", 
                           package = "panGenomeBreedr", 
                           mustWork = TRUE)

# Establish a virtual connection
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Fetch metadata for all sorghum accessions originating from Ethiopia
ethiopia_metadata_local <- fetch_accession_metadata(
  con = con_local,
  query_col = "countryorigin",
  query_value = "Ethiopia"
)
# print(ethiopia_metadata_local[1:3, 1:5])

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
