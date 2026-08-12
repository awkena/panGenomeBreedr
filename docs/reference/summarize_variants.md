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
# --- Online Mode ---
# Load the package
library(panGenomeBreedr)

# Get variant statistics from the online database, including annotation counts
online_stats <- summarize_variants(connect_db_mode = 'online', include_annotations = TRUE)
print(online_stats)
#>    chrom n_variants min_pos  max_pos n_unique_ids n_annotated
#> 1  Chr01    4179737     120 85112733      4179737     4179737
#> 2  Chr02    4197282      97 79114142      4197282     4197282
#> 3  Chr03    4067064     851 80862807      4067064     4067064
#> 4  Chr04    3921469     399 71213527      3921469     3921469
#> 5  Chr05    4503569     544 77058027      4503569     4503569
#> 6  Chr06    3362692     142 62713669      3362692     3362692
#> 7  Chr07    3709045     665 68910894      3709045     3709045
#> 8  Chr08    3540234     324 65779151      3540234     3540234
#> 9  Chr09    3265682   16102 63277331      3265682     3265682
#> 10 Chr10    3484137     424 62863324      3484137     3484137

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Get variant statistics across all chromosomes
local_stats <- summarize_variants(con = con_local, include_annotations = TRUE)
# print(local_stats)

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
