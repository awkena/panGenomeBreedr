# Fetch variants by allele frequency

This function queries variants within a specified genomic region and
filters them based on a minimum and maximum alternate allele frequency.
It can operate in either 'local' mode (querying a local DuckDB database)
or 'online' mode (fetching data from a remote API endpoint).

## Usage

``` r
fetch_variants_by_allele_frequency(
  con = NULL,
  min_af = 0,
  max_af = 1,
  chrom = NULL,
  start = NULL,
  end = NULL,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- min_af:

  Numeric. Minimum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `0`.

- max_af:

  Numeric. Maximum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `1`.

- chrom:

  A character value specifying the target chromosome name (e.g.,
  `"Chr05"`).

- start:

  Integer. Start coordinate for the target window region.

- end:

  Integer. End coordinate for the target window region.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame of variants passing the frequency thresholds within the
defined region, tracking their coordinates, calculated reference allele
frequencies (`ref_af`), and alternate allele frequencies (`alt_af`).

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# --- Online Mode ---
# Query for common variants (MAF >= 5%) in a specific region from the online database
online_common_vars <- fetch_variants_by_allele_frequency(
  min_af = 0.05,
  max_af = 0.95,
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)
print(head(online_common_vars))
#>             variant_id chrom      pos    ref_af     alt_af
#> 2 INDEL_Chr05_75104564 Chr05 75104564 0.8854415 0.11455847
#> 3 INDEL_Chr05_75104573 Chr05 75104573 0.8848449 0.11515513
#> 4 INDEL_Chr05_75104585 Chr05 75104585 0.8836516 0.11634845
#> 5 INDEL_Chr05_75104618 Chr05 75104618 0.8845465 0.11545346
#> 7 INDEL_Chr05_75104704 Chr05 75104704 0.9334726 0.06652745
#> 9 INDEL_Chr05_75104800 Chr05 75104800 0.8135442 0.18645585

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Query a specific locus window and filter for common variants (MAF >= 5%)
local_common_vars <- fetch_variants_by_allele_frequency(con = con_local,
                                min_af = 0.05,
                                max_af = 0.95,
                                chrom = "Chr05",
                                start = 75104537,
                                end = 75106403)
# print(head(local_common_vars))

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
