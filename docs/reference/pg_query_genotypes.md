# Query genotypes for specific variant IDs

This function connects to the public panGenomeBreedr API to retrieve
genomic data for a specific list of variant IDs. It expands the genotype
array into a wide format (samples as columns).

## Usage

``` r
pg_query_genotypes(
  variant_ids,
  variant_id_col = "variant_id",
  variants_table = "variants",
  genotypes_table = "genotypes",
  meta_data = NULL
)
```

## Arguments

- variant_ids:

  A character vector of variant IDs to retrieve.

- variant_id_col:

  Character. Name of the ID column. Default is "variant_id".

- variants_table:

  Character. Name of the metadata table. Default is "variants".

- genotypes_table:

  Character. Name of the genotype table. Default is "genotypes".

- meta_data:

  Character vector. Specific columns to include from the variants table
  (e.g., "chrom", "pos", "ref", "alt"). If `NULL`, retrieves all
  columns.

## Value

A data frame in wide format (variants x samples) with the requested
metadata columns.
