# List all tables in the SQLite database.

List all tables in the SQLite database.

## Usage

``` r
list_sqlite_tables(db_path)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

## Value

A character vector of names of SQLite databases

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# Define tempdir
path <- tempdir()

# Mini SQLite database
mini_db <-  system.file("extdata", "mini_sorghum_variant_vcf.db.gz",
                     package = "panGenomeBreedr",
                     mustWork = TRUE)

# Path to SQLite databases: INDEL and SNP
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')

# Unzip compressed mini database and save in tempdir
R.utils::gunzip(mini_db,
               destname = mini_db_path,
               remove = FALSE)

# List tables in the SQLite database
list_sqlite_tables(mini_db_path)
#> [1] "annotations" "genotypes"   "variants"   

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
