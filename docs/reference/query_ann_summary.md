# Query the annotations table within a specified genomic region and summarize the distribution of SnpEff annotations and impact categories by variant type.

Query the annotations table within a specified genomic region and
summarize the distribution of SnpEff annotations and impact categories
by variant type.

## Usage

``` r
query_ann_summary(
  db_path,
  chrom,
  start,
  end,
  annotations_table = "annotations",
  variants_table = "variants"
)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

- chrom:

  A character value specifying the chromosome name.

- start, end:

  A numeric value specifying the start and end coordinates for the
  candidate gene.

- annotations_table, variants_table:

  A character value indicating the table names for snpEff annotations
  and variants in the SQLite database.

## Value

A list with:

- `annotation_summary`

- `impact_summary`

- `variant_type_totals`

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

# Annotation and impact summary distribution within Chr05: 75104537-75106403
ann_summary <- query_ann_summary(db_path = mini_db_path,
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
