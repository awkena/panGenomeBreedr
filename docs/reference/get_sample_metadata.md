# Retrieve sample metadata from the SQLite database

This function connects to the local SQLite database and retrieves
accession-level metadata from the 'metadata' table. It supports optional
filtering by specific columns (e.g., country of origin, race, or
cluster) to easily subset populations for downstream analysis.

## Usage

``` r
get_sample_metadata(db_path, query_col = NULL, query_value = NULL)
```

## Arguments

- db_path:

  A character string specifying the file path to the SQLite database.

- query_col:

  Character. The metadata column to filter by (e.g., "countryorigin").
  If `NULL`, all records are returned.

- query_value:

  Character or Numeric. The specific value to match in `query_col`.

## Value

A data frame containing the sample metadata records, ordered by their
array index to seamlessly match genotype matrices.

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Define tempdir and setup mini SQLite database
path <- tempdir()
mini_db <- system.file(
  "extdata", "mini_sorghum_variant_vcf.db.gz", 
  package = "panGenomeBreedr", mustWork = TRUE
)
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)

# Fetch metadata for all accessions originating from Ethiopia
eth_samples <- get_sample_metadata(
  db_path = mini_db_path,
  query_col = "countryorigin",
  query_value = "Ethiopia"
)

# Clean tempdir
unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE, 
                  include.dirs = TRUE), recursive = TRUE, force = TRUE)
} # }
```
