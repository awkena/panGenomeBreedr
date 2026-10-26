# Get variant statistics from the database

Retrieves variant statistics, such as counts and genomic range, grouped
by chromosome. This function can operate in two modes: 'local' to query
a local database connection, or 'online' to fetch data from the API.

## Usage

``` r
summarize_variants(
  con = NULL,
  connect_db_mode = c("local", "online"),
  include_annotations = TRUE
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

- include_annotations:

  A logical value indicating whether to include statistics for the
  annotations table. Defaults to `TRUE`.

## Value

A data frame of variant statistics grouped by chromosome.

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

# Get variant statistics across all chromosomes
local_stats <- summarize_variants(con = con_local, include_annotations = TRUE)
print(local_stats)
#>   chrom n_variants  min_pos  max_pos n_unique_ids n_annotated
#> 1 Chr03         28 79037682 79039059           28          28
#> 2 Chr05        102 75104541 75106383          102         102

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To query the public online resource instead:
online_stats <- summarize_variants(connect_db_mode = 'online', include_annotations = TRUE)
print(online_stats)
} # }
```
