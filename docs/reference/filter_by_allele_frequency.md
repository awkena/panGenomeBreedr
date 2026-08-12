# Filter genotypes by allele frequency

Filter genotypes by allele frequency

## Usage

``` r
filter_by_allele_frequency(
  gt,
  variant_id_col = "variant_id",
  chrom_col = "chrom",
  pos_col = "pos",
  min_af = 0,
  max_af = 1
)
```

## Arguments

- gt:

  A data frame containing genotype calls queried from the genotypes
  table.

- variant_id_col:

  Character. The column name in `gt` matching unique variant
  identifiers. Defaults to `"variant_id"`.

- chrom_col:

  Character. The column name in `gt` matching chromosome names. Defaults
  to `"chrom"`.

- pos_col:

  Character. The column name in `gt` matching genomic positions.
  Defaults to `"pos"`.

- min_af:

  Numeric. Minimum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `0`.

- max_af:

  Numeric. Maximum alternate allele frequency threshold for filtering.
  Must fall between 0 and 1. Defaults to `1`.

## Value

A data frame of variants passing the frequency thresholds, tracking
their coordinates, calculated reference allele frequencies (`ref_af`),
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

# Query genotypes
gt_region <- fetch_table_region(con = con, chrom = "Chr05", start = 75104537, 
                      end = 75106403, table_name = "genotypes")

filtered_vars <- filter_by_allele_frequency(gt_region, variant_id_col = "variant_id", 
                              chrom_col = "chrom", pos_col = "pos", 
                              min_af = 0.05, max_af = 0.95)

# Disconnect at the end of your session
disconnect_local_db(con)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
