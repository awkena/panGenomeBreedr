# Get variants statistics stored in SQLite database

Get variants statistics stored in SQLite database

## Usage

``` r
variant_stats(db_path, include_annotations = TRUE)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

- include_annotations:

  A logical value indicating whether to include statistics for the
  annotations table.

## Value

A data frame of variant statistics in the SQLite databases

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

# Get variant statistics in the SQLite database
variant_stats(mini_db_path)
#>   chrom n_variants  min_pos  max_pos n_unique_ids n_annotated
#> 1 Chr05        216 75103571 75107320          216         216

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
