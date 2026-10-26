# Fetch Pangenome-Characterized Introgression Lines (PCIL) by Variant

This function identifies Pangenome-Characterized Introgression Lines
(PCIL) that are putative donors for a given variant of interest. It
operates by comparing the genotypes of recurrent parents against a pool
of potential donors to find families that carry the alternate allele.

## Usage

``` r
fetch_pcil_families_by_variant(
  con = NULL,
  selection,
  pcil_data,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database.

- selection:

  A character vector of variant IDs.

- pcil_data:

  A list object containing PCIL data, typically from
  \`fetch_pcil_data()\`. Must contain the \`pcil_metadata\` table.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A list containing summaries of PCIL families, genotypes, and selections.

- `pcil_family_summary`: A summary of PCIL families, including counts of
  families and lines per recurrent parent that are putative sources for
  the alternate allele.

- `geno_pi`: Genotype information for the parental lines.

- `pcil_summary`: Detailed information on the selected PCIL families.

Returns \`NULL\` if no matching PCIL families are found.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# 1. Connect to the package's bundled example database and load PCIL data
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)
con <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.
pcil_data <- fetch_pcil_data(con = con, connect_db_mode = "local")

# 2. Define variants of interest
selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")

# 3. Fetch PCIL families acting as putative donors for these variants
results <- fetch_pcil_families_by_variant(
  con = con,
  selection = selection,
  pcil_data = pcil_data,
  connect_db_mode = "local"
)

print(results$pcil_family_summary)
#>          clan families lines rp_genotype            selection recurrent_allele
#> 1     IRAT204        2    29   Reference INDEL_Chr03_79037889               TG
#> 2     IRAT204        2    29   Reference   SNP_Chr03_79037855                C
#> 3 Mota Maradi        1    18   Reference INDEL_Chr03_79037889               TG
#> 4 Mota Maradi        1    18   Reference   SNP_Chr03_79037855                C

disconnect_local_db(con)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
