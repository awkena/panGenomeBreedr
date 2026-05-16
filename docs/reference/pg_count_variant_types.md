# Count the number of variant types in the pangenome database.

This function connects to the public panGenomeBreedr API to perform a
server-side aggregation, counting the occurrences of different variant
types (e.g., SNP, INDEL) stored in the database.

## Usage

``` r
pg_count_variant_types(variants_table = "variants")
```

## Arguments

- variants_table:

  Character. The name of the table containing variant metadata. Defaults
  to "variants".

## Value

A data frame with two columns: `variant_type` and `n`.
