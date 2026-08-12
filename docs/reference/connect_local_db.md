# Connect to the local offline database

Establishes a virtual connection to the offline database engine using
DuckDB. It mounts the core pangenome-scale Parquet files .

## Usage

``` r
connect_local_db(folder_path, max_memory = "8GB", quiet = FALSE)
```

## Arguments

- folder_path:

  Character. The absolute or relative path to the directory containing
  the core database files (\`genotypes_partitioned/\`,
  \`variants.parquet\`, \`metadata.parquet\`, \`annotations.parquet\`)
  and optionally the \`pcil/\` extension directory.

- max_memory:

  Character. Maximum RAM DuckDB is allowed to utilize. Defaults to
  '"8GB"'. Adjust based on system constraints.

- quiet:

  Logical. If `TRUE`, suppresses console messages. Defaults to `FALSE`.

## Value

A DBI connection object (\`con\`) linked to the virtual database.

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)
my_db_folder <- system.file("extdata", "pangenome_scale_db", 
                           package = "panGenomeBreedr", 
                           mustWork = TRUE)


# Establish the connection
con <- connect_local_db(folder_path = my_db_folder, max_memory = "8GB")
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# }
```
