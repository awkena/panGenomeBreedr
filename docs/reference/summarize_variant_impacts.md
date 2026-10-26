# Get summary of database variant impacts by chromosome

This function summarizes the counts of different variant impact levels
(e.g., HIGH, MODERATE, LOW) for each chromosome and returns the result
in a wide-format data frame. It can operate in 'local' or 'online' mode.

## Usage

``` r
summarize_variant_impacts(con = NULL, connect_db_mode = c("local", "online"))
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A wide-format data frame where the first column is `chrom`, followed by
pivoted columns for each impact type (e.g., `impact_HIGH`,
`impact_MODERATE`), containing their respective mutation counts.

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

# Generate the wide-format impact profile matrix
local_impact <- summarize_variant_impacts(con = con_local)
print(local_impact)
#>   chrom impact_HIGH impact_LOW impact_MODERATE impact_MODIFIER
#> 1 Chr03           1         15               8              34
#> 2 Chr05           6         21              34             302

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To query the public online resource instead:
online_impact <- summarize_variant_impacts(connect_db_mode = 'online')
print(online_impact)
} # }
```
