# Summarize functional annotations and impacts for a genomic region

This function queries variants within a specific genomic range and
returns summaries of SnpEff annotations and impact levels,
cross-tabulated by variant type. It can operate in 'local' or 'online'
mode.

## Usage

``` r
summarize_annotations(
  con = NULL,
  chrom,
  start,
  end,
  annotations_table = "annotations",
  variants_table = "variants",
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- chrom:

  Character. Chromosome name (e.g., "Chr05").

- start:

  Numeric. Start coordinate of the region.

- end:

  Numeric. End coordinate of the region.

- annotations_table:

  Character. Name of the annotations table. Defaults to "annotations".

- variants_table:

  Character. Name of the variants table. Defaults to "variants".

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A list object containing three summarized summary data frames:

- `annotation_summary`: Counts of distinct functional mutations grouped
  by type.

- `impact_summary`: Counts of predicted functional impacts grouped by
  type.

- `variant_type_totals`: Entire marker breakdown counts across the
  defined coordinate window.

## Examples

``` r
# \donttest{
# --- Online Mode ---
# Load the package
library(panGenomeBreedr)

# Get annotation summary from the online database
online_summary <- summarize_annotations(
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)
print(online_summary)
#> $annotation_summary
#>                        annotation variant_type count
#> 1   conservative_inframe_deletion        INDEL     1
#> 2                  intron_variant        INDEL     2
#> 3             5_prime_UTR_variant          SNP     1
#> 4    disruptive_inframe_insertion        INDEL     1
#> 5  conservative_inframe_insertion        INDEL     1
#> 6           upstream_gene_variant          SNP   149
#> 7             3_prime_UTR_variant        INDEL    13
#> 8                missense_variant          SNP    29
#> 9           upstream_gene_variant        INDEL    53
#> 10             synonymous_variant          SNP    21
#> 11        downstream_gene_variant        INDEL    13
#> 12    disruptive_inframe_deletion        INDEL     2
#> 13             frameshift_variant        INDEL     6
#> 14            3_prime_UTR_variant          SNP    18
#> 15            5_prime_UTR_variant        INDEL     1
#> 16                 intron_variant          SNP     6
#> 17        downstream_gene_variant          SNP    46
#> 
#> $impact_summary
#>     impact variant_type count
#> 1 MODIFIER          SNP   220
#> 2 MODIFIER        INDEL    82
#> 3 MODERATE        INDEL     5
#> 4      LOW          SNP    21
#> 5 MODERATE          SNP    29
#> 6     HIGH        INDEL     6
#> 
#> $variant_type_totals
#>   variant_type total_variants
#> 1        INDEL             27
#> 2          SNP             75
#> 

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db", 
                           package = "panGenomeBreedr", 
                           mustWork = TRUE)

# Establish a virtual connection
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Extract functional distribution summaries across a specific gene locus window
local_summary <- summarize_annotations(
  con = con_local,
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)
# print(local_summary)

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
