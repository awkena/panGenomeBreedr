# Query database tables by genomic coordinates (online)

This function connects to the public panGenomeBreedr API to retrieve
data from the variants, annotations, or genotypes tables based on a
specific chromosome and genomic range.

## Usage

``` r
pg_query_db(
  table_name = c("variants", "annotations", "genotypes"),
  chrom = NULL,
  start = NULL,
  end = NULL,
  gene_name = NULL
)
```

## Arguments

- table_name:

  A character value specifying the target view to query. Must be one of
  `"variants"`, `"annotations"`, or `"genotypes"`.

- chrom:

  A character value specifying the target chromosome name (e.g.,
  `"Chr05"`).

- start:

  Integer. Start coordinate for the target window region.

- end:

  Integer. End coordinate for the target window region.

- gene_name:

  A character value indicating the specific Sobic gene model ID.
  Utilized explicitly when subsetting the `"annotations"` layout.

## Value

A data frame containing the queried genomic data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Query variants table for a specific genomic region
variants_data <- pg_query_db(
  table_name = "variants",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)
print(head(variants_data))

# Query annotations table for a specific gene
annotations_data <- pg_query_db(
  table_name = "annotations",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  gene_name = "Sobic.005G213600"
)
print(head(annotations_data))

# Query genotypes table for the same region
genotypes_data <- pg_query_db(table_name = "genotypes", chrom = "Chr05",
                              start = 75104537, end = 75106403)
print(head(genotypes_data))
} # }
```
