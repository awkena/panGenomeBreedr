# Interactive geographic exploration of sorghum accessions (local)

Generates an interactive web map via `leaflet` displaying the geographic
distribution of sorghum lines based on latitude and longitude
coordinates. Points are dynamically colored by a user-specified metadata
category and display scrollable, auto-formatted HTML tooltips containing
all available metadata attributes.

## Usage

``` r
query_map_accessions(metadata, color_by = "countryorigin")
```

## Arguments

- metadata:

  A data frame containing sample metadata tracking information. Must
  include explicit numeric `"lat"` and `"lon"` columns alongside the
  column specified in `color_by`.

- color_by:

  Character. The specific metadata column name used to group and assign
  discrete color palette attributes to markers. Defaults to
  `"countryorigin"`.

## Value

A `leaflet` map object (htmlwidget) displaying interactive geographic
coordinate plots.

## Details

The function automatically strips out tracking lines containing missing
spatial coordinates. Rather than hardcoding fixed tooltip categories, it
dynamically extracts all non-coordinate metadata attributes present in
the provided matrix, formatting them into clear, scrollable HTML popup
panels for every marker.

**Dependency Note:** To keep core package requirements lean, `leaflet`
and `tools` are configured as suggested package attachments. If they are
absent from the local execution workspace, the script will halt cleanly
and prompt appropriate installation commands.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Connect to the database pipeline
con <- connect_local_db(folder_path = "~/Desktop/curated_sorghum_variant_resource")

# Fetch passport metadata details for all accessions matching the resource
sample_metadata <- get_sample_metadata(con = con)

# Generate an interactive geographic distribution plot colored by
# K-Means genetic cluster assignments
query_map_accessions(sample_metadata, color_by = "kmeans_cluster")

# Safely close out the database connection
disconnect_local_db(con)
} # }
```
