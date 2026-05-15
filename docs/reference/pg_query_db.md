# Query any table in the pangenome database using chromosome and a genomic position range.

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

  Character. One of "variants", "annotations", or "genotypes".

- chrom:

  Character. The chromosome name (e.g., "Chr05").

- start, end:

  Numeric. The genomic start and end positions.

- gene_name:

  Character. Optional Sobic ID to filter annotations (e.g.,
  "Sobic.005G213600").

## Value

A data frame containing the queried genomic data.
