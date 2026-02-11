# Query genotypes for one or more variant IDs from a wide-format genotype table.

Query genotypes for one or more variant IDs from a wide-format genotype
table.

## Usage

``` r
query_genotypes(
  db_path,
  variant_ids,
  variant_id_col = "variant_id",
  variants_table = "variants",
  genotypes_table = "genotypes",
  meta_data = NULL
)
```

## Arguments

- db_path:

  The path to your SQLite or PostgreSQL database.

- variant_ids:

  A character vector of variant IDs to query.

- variant_id_col:

  A character value specifying the column name of variant IDs in the
  SQLite database.

- variants_table, genotypes_table, :

  A character value specifying the column names of variants and
  genotypes tables, respectively, in the SQLite database. if the
  \`meta_data = NULL\`, it returns all the metadata in the variants
  table.

- meta_data:

  A character vector of metadata to include with the genotype data.

## Value

A data.frame of genotype data in variants x samples format with
metadata.

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

# Extract all genotypes for all samples for any set of variants
geno_filter <- query_genotypes(db_path = mini_db_path,
                               variant_ids = c("INDEL_Chr05_75104541",
                                               "SNP_Chr05_75104557"),
                               meta_data = c('chrom', 'pos', 'ref', 'alt',
                                             'variant_type'))
# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)
# }
```
