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
# --- Online Mode ---
# Load the package
library(panGenomeBreedr)

# Get variant impact summary from the online database
online_impact <- summarize_variant_impacts(connect_db_mode = 'online')
print(online_impact)
#>    chrom impact_HIGH impact_LOW impact_MODERATE impact_MODIFIER
#> 1  Chr01        6965     100643           94242         7684401
#> 2  Chr02        6455      80519           84276         6920733
#> 3  Chr03        5733      78427           75675         6859442
#> 4  Chr04        5137      73419           71810         6478860
#> 5  Chr05        6372      71008           87243         6382307
#> 6  Chr06        4091      53603           55277         5198904
#> 7  Chr07        4094      47808           50546         5345799
#> 8  Chr08        4449      53821           60643         5045650
#> 9  Chr09        3628      49357           51485         5001667
#> 10 Chr10        4248      54241           54677         5361457

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Generate the wide-format impact profile matrix
local_impact <- summarize_variant_impacts(con = con_local)
# print(local_impact)

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
