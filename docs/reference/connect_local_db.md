# Connect to local offline database

Initializes a high-performance, offline DuckDB database engine and maps
local Parquet files to virtual tables for instant querying. This
function creates a persistent connection that should be passed to all
downstream query functions.

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

# Connect to the database at the start of your session
my_db_folder <- "~/Desktop/curated_sorghum_variant_resource"
con <- connect_local_db(folder_path = my_db_folder, max_memory = "8GB")

# Pass the connection to downstream functions
summary_df <- summarize_sqlite_tables(con)

# Disconnect at the end of your session
disconnect_local_db(con)
} # }
```
