# Fetch data from database tables by genomic coordinates

This function queries variant, annotation, genotype, or metadata tables
based on genomic coordinates (chromosome, start, end) or a gene name. It
can operate in either 'local' mode (connecting to a local DuckDB
database) or 'online' mode (fetching data from a remote API endpoint).

## Usage

``` r
fetch_table_region(
  con = NULL,
  table_name = c("variants", "annotations", "genotypes"),
  chrom,
  start = NULL,
  end = NULL,
  gene_name = NULL,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- table_name:

  A character value specifying the target view to query. Must be one of
  \`"variants"\`, \`"annotations"\`, \`"genotypes"\`, or \`"metadata"\`.

- chrom:

  A character value specifying the target chromosome name (e.g.,
  \`"Chr05"\`).

- start:

  Integer. Optional start coordinate for the target window region.

- end:

  Integer. Optional end coordinate for the target window region.

- gene_name:

  A character value indicating the specific Sobic gene model ID.
  Utilized explicitly when subsetting the \`"annotations"\` table.

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame containing the targeted genomic records. For the
\`"genotypes"\` table, it returns individual samples unpacked directly
into intuitive wide columns, along with calculated allele metrics.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# --- Online Mode ---
# Query variants table for a specific genomic region from the online database
online_variants_data <- fetch_table_region(
  table_name = "variants",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)
print(head(online_variants_data))
#>           variant_id chrom      pos ref alt qual filter variant_type
#> 1 SNP_Chr05_75104557 Chr05 75104557   C   T    .   PASS          SNP
#> 2 SNP_Chr05_75104560 Chr05 75104560   C   T    .   PASS          SNP
#> 3 SNP_Chr05_75104568 Chr05 75104568   G   T    .   PASS          SNP
#> 4 SNP_Chr05_75104569 Chr05 75104569   C   A    .   PASS          SNP
#> 5 SNP_Chr05_75104574 Chr05 75104574   C   T    .   PASS          SNP
#> 6 SNP_Chr05_75104591 Chr05 75104591   A   G    .   PASS          SNP

# Query genotypes table for the same region from the online database
online_genotypes_data <- fetch_table_region(
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)
print(online_genotypes_data[1:6, 1:12])
#>           variant_id chrom      pos ref alt variant_type major_allele
#> 1 SNP_Chr05_75104557 Chr05 75104557   C   T          SNP            C
#> 2 SNP_Chr05_75104560 Chr05 75104560   C   T          SNP            C
#> 3 SNP_Chr05_75104568 Chr05 75104568   G   T          SNP            G
#> 4 SNP_Chr05_75104569 Chr05 75104569   C   A          SNP            C
#> 5 SNP_Chr05_75104574 Chr05 75104574   C   T          SNP            C
#> 6 SNP_Chr05_75104591 Chr05 75104591   A   G          SNP            A
#>   minor_allele major_allele_freq minor_allele_freq IDMM ISGC
#> 1            T           0.89051           0.10949  0|0  0|0
#> 2            T           0.88962           0.11038  0|0  0|0
#> 3            T           0.99791           0.00209  0|0  0|0
#> 4            A           0.99791           0.00209  0|0  0|0
#> 5            T           0.99195           0.00805  0|0  0|0
#> 6            G           0.96420           0.03580  0|0  0|0

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Extract functional annotations inside a candidate locus region from local DB
local_annota_region <- fetch_table_region(
  con = con_local,
  table_name = "annotations",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  gene_name = "Sobic.005G213600",
  connect_db_mode = 'local'
)
# print(head(local_annota_region))

# Extract matrix genotypes within the exact same coordinates window from local DB
local_gt_region <- fetch_table_region(
  con = con_local,
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'local'
)
# print(local_gt_region[1:6, 1:12])

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
