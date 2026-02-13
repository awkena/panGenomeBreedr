# Get variants statistics stored in SQLite database based on mutation impact.

Get variants statistics stored in SQLite database based on mutation
impact.

## Usage

``` r
variant_impact_summary(db_path)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

## Value

A data frame of variant impact summary or counts.

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

# Get variant impact summary in the SQLite database
variant_impact_summary(mini_db_path)
#>   chrom impact_HIGH impact_LOW impact_MODERATE impact_MODIFIER
#> 1 Chr05           6         21              34             734

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
