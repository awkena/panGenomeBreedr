# Check the column names and types for any table in a SQLite database.

Check the column names and types for any table in a SQLite database.

## Usage

``` r
list_table_columns(
  db_path,
  table_name = c("variants", "annotations", "genotypes")
)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

- table_name:

  A character value specifying the name of any table in a SQLite
  database.

## Value

A data frame of table column information in the SQLite database.

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

# List column names and type of variable in the variants table in the database
list_table_columns(mini_db_path, table_name = "variants")
#>   cid         name    type notnull dflt_value pk
#> 1   0   variant_id    TEXT       0         NA  0
#> 2   1        chrom    TEXT       0         NA  0
#> 3   2          pos INTEGER       0         NA  0
#> 4   3          ref    TEXT       0         NA  0
#> 5   4          alt    TEXT       0         NA  0
#> 6   5         qual    TEXT       0         NA  0
#> 7   6       filter    TEXT       0         NA  0
#> 8   7 variant_type    TEXT       0         NA  0

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
