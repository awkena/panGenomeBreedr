# Summarize functional annotations and impacts for a genomic region (online)

This function connects to the public panGenomeBreedr API to query
variants within a specific genomic range and returns summaries of SnpEff
annotations and impact levels, cross-tabulated by variant type.

## Usage

``` r
pg_query_ann_summary(
  chrom,
  start,
  end,
  annotations_table = "annotations",
  variants_table = "variants"
)
```

## Arguments

- chrom:

  Character. Chromosome name (e.g., "Chr05").

- start:

  Numeric. Start coordinate of the region.

- end:

  Numeric. End coordinate of the region.

- annotations_table:

  Character. Name of the annotations table. Defaults to "annotations".

- variants_table:

  Character. Name of the variants table. Defaults to "variants".

## Value

A list containing three data frames: `annotation_summary`,
`impact_summary`, and `variant_type_totals`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Extract functional distribution summaries across a specific gene locus window
query_ann_summary_result <- pg_query_ann_summary(chrom = "Chr05",
                                              start = 75104537,
                                              end = 75106403)
print(query_ann_summary_result)

} # }
```
