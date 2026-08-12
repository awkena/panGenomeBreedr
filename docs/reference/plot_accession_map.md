# Interactive geographic exploration of the samples

Interactive geographic exploration of the samples

## Usage

``` r
plot_accession_map(metadata, color_by = "countryorigin")
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
# \donttest{
# Load the package
library(panGenomeBreedr)

# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db", 
                           package = "panGenomeBreedr", 
                           mustWork = TRUE)

# Establish a virtual connection to the offline database engine
con <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database!  Pangenome-scale database  mounted safely. No folder named pcil.

# Fetch passport metadata details for all accessions matching the resource
sample_metadata <- fetch_accession_metadata(con = con)

# Generate an interactive geographic distribution plot colored by countryorigin
plot_accession_map_result <- plot_accession_map(sample_metadata, color_by = "countryorigin")
print(plot_accession_map_result)

# Disconnect at the end of your session
disconnect_local_db(con)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
