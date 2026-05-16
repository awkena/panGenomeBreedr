# Retrieve sample metadata from the pangenome database.

This function connects to the public panGenomeBreedr API to fetch
accession-level metadata, such as origin, race, and classification from
the 'metadata' table. It supports optional filtering by specific columns
to easily subset populations for downstream analysis.

## Usage

``` r
pg_get_sample_metadata(query_col = NULL, query_value = NULL)
```

## Arguments

- query_col:

  Character. The metadata column to filter by (e.g., "countryorigin").
  If `NULL`, all records are returned.

- query_value:

  Character or Numeric. The specific value to match in `query_col`.

## Value

A data frame containing the sample metadata records.
