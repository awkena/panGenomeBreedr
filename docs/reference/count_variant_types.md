# Count the number of variant types in the SQLite database.

Count the number of variant types in the SQLite database.

## Usage

``` r
count_variant_types(db_path)
```

## Arguments

- db_path:

  The path to your SQLite or PostgreSQL database.

## Value

A data.frame of counts for variant types.

## Examples

``` r
# \donttest{
# example code
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

# Count the number of variant types in database
variant_type_count <- count_variant_types(mini_db_path)

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
