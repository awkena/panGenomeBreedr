# Interactive Geographic Exploration of Sorghum Accessions

This function generates a high-performance interactive map showing the
geographic distribution of sorghum lines. It dynamically generates rich,
scrollable popups containing all available metadata for each accession.

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

## Details

The function automatically filters out records with missing geographic
coordinates. Instead of hardcoding specific tooltip values, it
dynamically reads all non-coordinate columns from the provided
`metadata` data frame and formats them into a scrollable HTML popup for
each point.

**Dependency Note:** To keep the core package lightweight, the `leaflet`
and `tools` packages are listed as "Suggested" dependencies. If they are
not currently installed on your system, the function will gracefully
stop and prompt you to install them before plotting.

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Fetch sample metadata from the database
meta <- pg_get_sample_metadata()

# Explore the geographic distribution colored by genetic cluster
pg_map_accessions(meta, color_by = "kmeans_cluster")
} # }
```
