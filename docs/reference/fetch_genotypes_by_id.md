# Fetch genotypes for specific variant IDs

This function queries the genotypes table for a specific list of variant
IDs. It can operate in 'local' mode (querying a local DuckDB database)
or 'online' mode (fetching data from a remote API endpoint).

## Usage

``` r
fetch_genotypes_by_id(
  con = NULL,
  variant_ids,
  chrom = NULL,
  variant_id_col = "variant_id",
  variants_table = "variants",
  genotypes_table = "genotypes",
  meta_data = NULL,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- variant_ids:

  A character vector of unique variant identifiers to look up.

- chrom:

  An optional character vector of chromosome names to filter the search.
  If \`NULL\` (default), the function will attempt to automatically
  detect chromosome names from the \`variant_ids\` to optimize the
  query.

- variant_id_col:

  A character value specifying the primary key column name for variant
  IDs. Defaults to \`"variant_id"\`.

- variants_table:

  A character value specifying the name of the variants reference table.
  Defaults to \`"variants"\`.

- genotypes_table:

  A character value specifying the name of the genotype call data table.
  Defaults to \`"genotypes"\`.

- meta_data:

  A character vector of metadata columns to return. This can include any
  field from the \`variants\` table (e.g., \`"chrom"\`, \`"pos"\`,
  \`"ref"\`, \`"alt"\`) as well as dynamically calculated metrics:
  \`"major_allele"\`, \`"minor_allele"\`, \`"major_allele_freq"\`, and
  \`"minor_allele_freq"\`. If \`NULL\` (default), all available database
  fields and metrics are returned.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame containing the requested metadata and genotype calls in a
wide variants-by-samples matrix.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# Define a list of target variant IDs
target_markers <- c("INDEL_Chr05_75104541", "SNP_Chr05_75104557")

# --- Online Mode ---.
# Query genotypes from the online database
online_genotypes <- fetch_genotypes_by_id(
  variant_ids = target_markers, 
  meta_data = c("chrom", "pos", "ref", "alt", "minor_allele_freq"),
  connect_db_mode = 'online'
)
print(head(online_genotypes[,1:9]))
#>             variant_id chrom      pos ref  alt minor_allele_freq chrom_1
#> 1   SNP_Chr05_75104557 Chr05 75104557   C    T           0.10949   Chr05
#> 2 INDEL_Chr05_75104541 Chr05 75104541   T TGAC           0.00119   Chr05
#>      pos_1 IDMM
#> 1 75104557  0|0
#> 2 75104541  0|0

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Query genotypes from the local database
local_genotypes <- fetch_genotypes_by_id(
  con = con_local,
  variant_ids = target_markers,
  meta_data = c('chrom', 'minor_allele_freq')
)
# print(local_genotypes[,1:6])

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
