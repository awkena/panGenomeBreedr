# Interactive geographic exploration of the samples (local)

Interactive geographic exploration of the samples (local)

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Define the path to the curated_sorghum_variant_resource folder
my_db_folder <- "~/Desktop/curated_sorghum_variant_resource"

# Establish the connection
con <- connect_local_db(folder_path = my_db_folder)

# Fetch passport metadata details for all accessions matching the resource
sample_metadata <- get_sample_metadata(con = con)

# Generate an interactive geographic distribution plot colored by countryorigin
query_map_accessions_result <- query_map_accessions(sample_metadata, color_by = "countryorigin")
print(query_map_accessions_result)

# Disconnect at the end of your session
disconnect_local_db(con)
} # }
```
