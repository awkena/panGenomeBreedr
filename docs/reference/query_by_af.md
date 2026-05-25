# Extract variants based on minimum and maximum allele frequencies within a defined region in a SQLite database.

Extract variants based on minimum and maximum allele frequencies within
a defined region in a SQLite database.

## Usage

``` r
query_by_af(
  db_path,
  min_af = 0,
  max_af = 1,
  chrom = NULL,
  start = NULL,
  end = NULL
)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

- min_af, max_af:

  A numeric value indicating the minimum and maximum allele frequency
  thresholds for filtering variants.

- chrom:

  A character value specifying the chromosome name.

- start, end:

  A numeric value specifying the start and end coordinates for the
  candidate gene.

## Value

A data frame of filtered variants and their allele frequencies.

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

# Query database to filter variants in a region based on alt allele frequency
filter_af <- query_by_af(db_path = mini_db_path,
                         min_af = 0.01,
                         max_af = 0.99,
                         chrom = "Chr05",
                         start = 75104537,
                         end = 75106403)

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
