# Query any table in your SQLite database using chromosome and a genomic position range.

Query any table in your SQLite database using chromosome and a genomic
position range.

## Usage

``` r
query_db(
  db_path,
  table_name = c("variants", "annotations", "genotypes"),
  chrom,
  start = NULL,
  end = NULL,
  gene_name = NULL
)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

- table_name:

  A character value specifying the name of any any table in a SQLite
  database.

- chrom:

  A character value specifying the chromosome name.

- start, end:

  A numeric value specifying the start and end coordinates for the
  candidate gene.

- gene_name:

  A character value indicating the Sobic ID of candidate gene, if
  available. It used if querying the \`annotations\` table in the
  database.

## Value

A data frame of the table queried in the SQLite database.

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

# Extract snpEff annotation for variants within the candidate gene Sobic.005G213600
annota_region <- query_db(db_path = mini_db_path,
                           chrom = "Chr05",
                           start = 75104537,
                           end = 75106403,
                           table_name = "annotations",
                           gene_name = "Sobic.005G213600")

# Extract VCF genotypes for variants within the region: 'Chr05:75104537-75106403'
gt_region <- query_db(db_path = mini_db_path,
                      chrom = "Chr05",
                      start = 75104537,
                      end = 75106403,
                      table_name = "genotypes")

# Clean tempdir
contents <- list.files(tempdir(),
                             full.names = TRUE,
                             recursive = TRUE,
                             all.files = TRUE,
                             include.dirs = TRUE)
unlink(contents, recursive = TRUE, force = TRUE)

# }
```
