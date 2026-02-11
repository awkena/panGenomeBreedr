# Extract variants from annotation table based on impact type: LOW, MODERATE, HIGH, MODIFIER.

Extract variants from annotation table based on impact type: LOW,
MODERATE, HIGH, MODIFIER.

## Usage

``` r
query_by_impact(
  db_path,
  impact_level = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
  chrom = NULL,
  start = NULL,
  end = NULL
)
```

## Arguments

- db_path:

  A character value indicating the path to the SQLite database.

- impact_level:

  A character value specifying the vriant impact type following snpEff
  annotation conventions: \`HIGH\`, \`MODERATE\`, \`LOW\`, \`MODIFIER\`

- chrom:

  A character value specifying the chromosome name.

- start, end:

  A numeric value specifying the start and end coordinates for the
  candidate gene.

## Value

A data frame of variants showing the impact type defined in the SQLite
database.

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

# Extract low impact variant for a region or gene
high_variants <- query_by_impact(db_path = mini_db_path,
                               impact_level = 'high',
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
