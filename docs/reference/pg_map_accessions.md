# Interactive geographic exploration of the samples (online)

Interactive geographic exploration of the samples (online)

## Usage

``` r
pg_map_accessions(metadata, color_by = "countryorigin")
```

## Arguments

- metadata:

  A data frame containing sample metadata. Must include 'lat' and 'lon'
  columns, along with the specified coloring column.

- color_by:

  Character. The metadata column to use for point coloration. Defaults
  to "countryorigin".

## Value

A `leaflet` map object (htmlwidget) representing the interactive map.

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Fetch sample metadata from the database
meta <- pg_get_sample_metadata()

# Explore the geographic distribution colored by genetic cluster
pg_map_accessions(meta, color_by = "countryorigin")
} # }
```
