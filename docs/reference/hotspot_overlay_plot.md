# Generate a combined gene model and variant hotspot overlay plot (local)

Generates a vertically stacked visualization aligning a structural gene
model with its corresponding regional variant hotspots.

## Usage

``` r
hotspot_overlay_plot(
  gene_name,
  gff_path,
  connect_db_mode = c("online", "local"),
  local_db_path = NULL,
  selected_variants = NULL,
  text_sz = 2.5
)
```

## Arguments

- gene_name:

  A character string indicating the Sobic ID of the candidate gene.

- gff_path:

  A character string specifying the path to the GFF3 file.

- connect_db_mode:

  A character string specifying the database connection mode: either
  "online" or "local".

- local_db_path:

  A character string specifying the local folder path to the database if
  running in "local" mode.

- selected_variants:

  An optional character vector of specific \`variant_id\`s to highlight
  on the plot with vertical lines, enlarged points, and labels.

- text_sz:

  A numeric value for specifying text size for selected variants.

## Value

A `patchwork` object containing the combined, aligned plots.

## Details

This wrapper function integrates genomic annotations and variant
genotypes. It natively queries the pangenome-scale variant relational
database. Whether operating in 'online' mode or 'local' mode, it
leverages the underlying Parquet and DuckDB architectures optimized for
sorghum to efficiently retrieve variants and their snpEff annotations
before plotting.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)
gff_path <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
hotspot_overlay_plot(gene_name = "Sobic.005G213600",
                     gff_path = gff_path,
                     connect_db_mode = "online",
                     selected_variants = c("INDEL_Chr05_75104881", "INDEL_Chr05_75105587"))

# }
```
