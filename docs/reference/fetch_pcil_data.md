# Fetch All PCIL Data Tables

Fetches all tables related to Pangenome-Characterized Introgression
Lines (PCIL) and returns them as a single list object. This function
serves as a "load all" utility for PCIL data, supporting both local and
online database modes.

## Usage

``` r
fetch_pcil_data(con = NULL, connect_db_mode = c("local", "online"))
```

## Arguments

- con:

  A DBI connection object to the local database. Required when
  \`connect_db_mode = 'local'\`.

- connect_db_mode:

  Character string specifying the connection mode: \`'local'\` or
  \`'online'\`.

## Value

A list containing the following data frames:

- `pcil_metadata`: Metadata for PCIL families.

- `pcil_gene_regions`: Genomic regions of genes within the PCILs.

- `pcil_introgressions`: Introgression block data.

- `pcil_genomewide_introgressions`: Genome-wide introgression summary
  statistics.

- `pcil_inbreeding_coefficient`: Sample metadata, often used for
  inbreeding coefficients.

- `pcil_IBS_dis`: Identity-by-state (IBS) distance matrix.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# --- Online Mode ---
# Fetch all PCIL data tables from the remote database
pcil_data <- fetch_pcil_data(connect_db_mode = 'online')
# }

if (FALSE) { # \dontrun{
# --- Offline Mode ---
# Requires a local database folder that contains a "pcil" subfolder;
# the sample dataset shipped with this package does not include one.
con <- connect_local_db(folder_path = "path/to/your/local_db")
pcil_data_local <- fetch_pcil_data(con = con, connect_db_mode = "local")
disconnect_local_db(con)
} # }
```
