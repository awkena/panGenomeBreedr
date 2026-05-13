# Query genotypes filtered by sample metadata attributes in SQLite

This function retrieves genotypes for a specific genomic region, but
restricts the output to a subset of samples defined by metadata
attributes (e.g., extracting variants only for specific countries,
populations, or phenotypic clusters).

## Usage

``` r
query_by_metadata(db_path, chrom, start, end, meta_col, meta_value)
```

## Arguments

- db_path:

  A character string specifying the file path to the SQLite database.

- chrom:

  Character. Chromosome name (e.g., "Chr05").

- start, end:

  Numeric. The genomic start and end coordinates.

- meta_col:

  Character. Metadata column to filter samples by (e.g.,
  "countryorigin").

- meta_value:

  Character. Value to match in `meta_col` (e.g., "Mali").

## Value

A wide-format data frame (variants x filtered samples).

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Define tempdir and setup mini SQLite database
path <- tempdir()
mini_db <- system.file(
  "extdata", "mini_sorghum_variant_vcf.db.gz", 
  package = "panGenomeBreedr", mustWork = TRUE
)
mini_db_path <- file.path(path, 'mini_sorghum_variant_vcf.db')
R.utils::gunzip(mini_db, destname = mini_db_path, remove = FALSE)

# Get genotypes for a specific gene region, but only for samples from Mali
mali_genotypes <- query_by_metadata(
  db_path = mini_db_path,
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  meta_col = "countryorigin",
  meta_value = "Mali"
)

# Clean tempdir
unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE, 
                  include.dirs = TRUE), recursive = TRUE, force = TRUE)
} # }
```
