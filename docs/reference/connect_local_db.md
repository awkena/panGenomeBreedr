# Connect to the local offline database

Connect to the local offline database

## Usage

``` r
connect_local_db(folder_path, max_memory = "8GB", quiet = FALSE)
```

## Arguments

- folder_path:

  Character. The absolute or relative path to the directory containing
  the core Parquet files (\`genotypes.parquet\`, \`variants.parquet\`,
  \`metadata.parquet\`, \`annotations.parquet\`).

- max_memory:

  Character. Maximum RAM DuckDB is allowed to utilize. Defaults to
  '"8GB"'. Adjust based on system constraints.

- quiet:

  Logical. If `TRUE`, suppresses console messages. Defaults to `FALSE`.

## Value

A DBI connection object (\`con\`) linked to the virtual database.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Define the path to the curated_sorghum_variant_resource folder
my_db_folder <- "~/Desktop/curated_sorghum_variant_resource"

# Establish the connection
con <- connect_local_db(folder_path = my_db_folder, max_memory = "8GB")

} # }
```
