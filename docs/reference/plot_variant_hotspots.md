# Plot variant hotspots

Plot variant hotspots

## Usage

``` r
plot_variant_hotspots(
  variant_table,
  annotation_table = NULL,
  region_start = NULL,
  region_end = NULL
)
```

## Arguments

- variant_table:

  Data frame containing basic variants

- annotation_table:

  Data frame containing annotations

- region_start:

  Numeric start position

- region_end:

  Numeric end position

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# Fetch variant and annotation data for a specific region from the online database
variants <- fetch_table_region(
  table_name = "variants",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)

annotations <- fetch_table_region(
  table_name = "annotations",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403,
  connect_db_mode = 'online'
)

# Generate the variant hotspot plot
hotspot_plot <- plot_variant_hotspots(
  variant_table = variants,
  annotation_table = annotations
)

# Display the plot
print(hotspot_plot)

# }
```
