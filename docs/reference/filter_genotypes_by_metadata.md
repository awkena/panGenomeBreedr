# Filter a genotype matrix by sample metadata

This function filters a genotype matrix based on sample metadata
criteria, allowing for selection of specific sub-populations. It can
operate in 'local' mode (querying a local DuckDB database for metadata)
or 'online' mode (fetching metadata from a remote API endpoint).

## Usage

``` r
filter_genotypes_by_metadata(
  con = NULL,
  genotype_matrix,
  genotype_start_col,
  filters,
  connect_db_mode = c("local", "online")
)
```

## Arguments

- con:

  A DBI connection object to the local database. Required only when
  \`connect_db_mode\` is 'local'. Defaults to \`NULL\`.

- genotype_matrix:

  A data frame or matrix where rows are variants and columns include
  variant metadata followed by sample genotypes.

- genotype_start_col:

  An integer specifying the column index where sample genotype calls
  begins.

- filters:

  A named list of metadata criteria to filter samples by. Names must
  match metadata schema columns, and values are the allowed entries.
  Example:
  `list(population = "Gates", countryorigin = c("Ethiopia", "Ghana", "Togo"))`

- connect_db_mode:

  A character string specifying the connection mode. Can be either
  \`'local'\` (default) or \`'online'\`.

## Value

A data frame containing the filtered genotype matrix with recalculated
allele metrics for the sub-population. Returns an empty data frame if no
samples match all combined criteria.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# Define filtering criteria
my_filters <- list(
  population = "Gates",
  countryorigin = c("Ethiopia", "Ghana", "Togo")
)

# --- Online Mode ---
# Extract the full genotype matrix first for a genomic region from the online database
online_genotype_data_region <- fetch_table_region(
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)

# Get filtered genotypes matrix from online data
online_filtered_genotypes <- filter_genotypes_by_metadata(
  genotype_matrix = online_genotype_data_region,
  genotype_start_col = 11, # Assuming 11 is the correct start column for this data
  filters = my_filters,
  connect_db_mode = 'online'
)
print(head(online_filtered_genotypes))
#>             variant_id chrom      pos variant_type ref  alt major_allele
#> 1 INDEL_Chr05_75104541 Chr05 75104541        INDEL   T TGAC            T
#> 2 INDEL_Chr05_75104564 Chr05 75104564        INDEL   C   CA            C
#> 3 INDEL_Chr05_75104573 Chr05 75104573        INDEL  TC    T           TC
#> 4 INDEL_Chr05_75104585 Chr05 75104585        INDEL   A  AAT            A
#> 5 INDEL_Chr05_75104618 Chr05 75104618        INDEL   G   GC            G
#> 6 INDEL_Chr05_75104703 Chr05 75104703        INDEL   C  CTG            C
#>   minor_allele major_allele_freq minor_allele_freq JHGC ISGD IQSD ISHD JBSS
#> 1         TGAC            1.0000            0.0000  0|0  0|0  0|0  0|0  0|0
#> 2           CA            0.7234            0.2766  0|0  0|0  0|0  1|1  0|0
#> 3            T            0.7234            0.2766  0|0  0|0  0|0  1|1  0|0
#> 4          AAT            0.7234            0.2766  0|0  0|0  0|0  1|1  0|0
#> 5           GC            0.7234            0.2766  0|0  0|0  0|0  1|1  0|0
#> 6          CTG            1.0000            0.0000  0|0  0|0  0|0  0|0  0|0
#>   JBST IZJQ IZJR IZJS IZJT IZJU IZJW IZJX IZJY IZJZ IZKA INUQ INUT IQSU IQSW
#> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
#> 2  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1
#> 3  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1
#> 4  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1
#> 5  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  1|1  1|1
#> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
#>   IQST INUX INXE INXF JHGT IXPJ IXNM INXH INXI IXML IXPK INXJ IZKP IXNN IXNP
#> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
#> 2  0|0  1|1  0|0  0|0  1|1  1|1  0|0  0|0  0|0  1|1  1|1  0|0  0|0  1|1  1|1
#> 3  0|0  1|1  0|0  0|0  1|1  1|1  0|0  0|0  0|0  1|1  1|1  0|0  0|0  1|1  1|1
#> 4  0|0  1|1  0|0  0|0  1|1  1|1  0|0  0|0  0|0  1|1  1|1  0|0  0|0  1|1  1|1
#> 5  0|0  1|1  0|0  0|0  1|1  1|1  0|0  0|0  0|0  1|1  1|1  0|0  0|0  1|1  1|1
#> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
#>   IPCI IPCJ ISHX ISIF ISIN ISIX ISJF ISHY JBPS JHHC IXNH JHHD
#> 1  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0
#> 2  0|0  0|0  0|0  0|0  1|1  0|0  1|1  1|1  0|0  0|0  0|0  0|0
#> 3  0|0  0|0  0|0  0|0  1|1  0|0  1|1  1|1  0|0  0|0  0|0  0|0
#> 4  0|0  0|0  0|0  0|0  1|1  0|0  1|1  1|1  0|0  0|0  0|0  0|0
#> 5  0|0  0|0  0|0  0|0  1|1  0|0  1|1  1|1  0|0  0|0  0|0  0|0
#> 6  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0  0|0

# --- Offline Mode ---
# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db", 
                           package = "panGenomeBreedr", 
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Extract the full genotype matrix first for the genomic region from local DB
local_genotype_data_region <- fetch_table_region(
  con = con_local,
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)

# Get filtered genotypes matrix from local data
local_filtered_genotypes <- filter_genotypes_by_metadata(
  con = con_local,
  genotype_matrix = local_genotype_data_region,
  genotype_start_col = 11,
  filters = my_filters
)
# print(head(local_filtered_genotypes))

# Disconnect at the end of your session
disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
