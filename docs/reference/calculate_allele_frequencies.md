# Compute allele frequencies for a genotype matrix

Compute allele frequencies for a genotype matrix

## Usage

``` r
calculate_allele_frequencies(
  gt,
  variant_id_col = "variant_id",
  chrom_col = NULL,
  pos_col = NULL
)
```

## Arguments

- gt:

  A matrix or data frame containing variants (rows) across samples
  (columns), including required variant identification and optional
  genomic coordinate metadata columns.

- variant_id_col:

  Character. The column name in `gt` matching unique variant
  identifiers. Defaults to `"variant_id"`.

- chrom_col:

  Character. Optional column name in `gt` matching chromosome names.
  Defaults to `NULL`.

- pos_col:

  Character. Optional column name in `gt` matching genomic positions.
  Defaults to `NULL`.

## Value

A data frame containing unique variant identifiers, optional coordinate
tracking columns, calculated reference allele frequencies (`ref_af`),
and alternate allele frequencies (`alt_af`).

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db", 
                           package = "panGenomeBreedr", 
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Pull wide genotype records for a locus window range
gt_region <- fetch_table_region(con = con, chrom = "Chr05", start = 75104537, 
                      end = 75106403, table_name = "genotypes")

# Compute allele frequencies across all samples in the dataset.
af_metrics <- calculate_allele_frequencies(gt_region, variant_id_col = "variant_id", 
                      chrom_col = "chrom", pos_col = "pos")

# Disconnect at the end of your session
disconnect_local_db(con)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
