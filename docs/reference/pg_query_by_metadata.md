# Query genotypes filtered by sample metadata attributes in the pangenome database.

This function connects to the public panGenomeBreedr API to retrieve
genotypes for a specific genomic region, but restricts the output to a
subset of samples defined by metadata attributes (e.g., extracting
variants only for specific countries, populations, or phenotypic
clusters).

## Usage

``` r
pg_query_by_metadata(chrom, start, end, meta_col, meta_value)
```

## Arguments

- chrom:

  Character. Chromosome name.

- start:

  Numeric. Genomic start range.

- end:

  Numeric. Genomic end range.

- meta_col:

  Character. Metadata column to filter samples by (e.g.,
  "countryorigin").

- meta_value:

  Character. Value to match in `meta_col` (e.g., "Ghana").

## Value

A wide-format data frame (variants x filtered samples).

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Get genotypes for a gene, but only for samples from Ethiopia
eth_genotypes <- pg_query_by_metadata(
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  meta_col = "countryorigin",
  meta_value = "Ethiopia"
)
} # }
```
