# Generate a combined gene model and variant hotspot overlay plot (local)

Generates a vertically stacked visualization aligning a structural gene
model with its corresponding regional variant hotspots.

## Usage

``` r
hotspot_overlay_plot(
  gene_name,
  gff_path,
  annotations_df,
  genotypes_df,
  selected_variants = NULL,
  text_sz = 2.5
)
```

## Arguments

- gene_name:

  A character string indicating the Sobic ID of the candidate gene.

- gff_path:

  A character string specifying the path to the GFF3 file.

- annotations_df:

  A data frame containing variant annotations for the region.

- genotypes_df:

  A data frame containing variant genotypes for the region.

- selected_variants:

  An optional character vector of specific \`variant_id\`s to highlight
  on the plot with vertical lines, enlarged points, and labels.

- text_sz:

  A numeric value for specifying text size for selected variants.

## Value

A `patchwork` object containing the combined, aligned plots.

## Details

This wrapper function integrates genomic annotations and variant
genotypes.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# 1. Define parameters using the package's bundled example database and
# a small GFF3 slice for the same gene
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)
gff_path <- system.file("extdata", "pangenome_scale_db", "gene_models.gff3.gz",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.
gene <- "Sobic.005G213600"

# Fetch data for the gene region
ann_df <- fetch_table_region(
  con = con_local, table_name = "annotations",
  chrom = "Chr05", start = 75104537, end = 75106403
)
geno_df <- fetch_table_region(
  con = con_local, table_name = "genotypes",
  chrom = "Chr05", start = 75104537, end = 75106403
)

#  Generate and display the plot
if (nrow(ann_df) > 0 && nrow(geno_df) > 0) {
  hotspot_overlay_plot(
    gene_name = gene, gff_path = gff_path,
    annotations_df = ann_df, genotypes_df = geno_df
  )
}


disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To use the full, public reference GFF3 and query the online resource instead:
gff_path <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
gene <- "Sobic.005G213600"
ann_df_online <- fetch_table_region(
  table_name = "annotations",
  chrom = "Chr05", start = 75104537, end = 75106403,
  connect_db_mode = 'online'
)
geno_df_online <- fetch_table_region(
  table_name = "genotypes",
  chrom = "Chr05", start = 75104537, end = 75106403,
  connect_db_mode = 'online'
)
if (nrow(ann_df_online) > 0 && nrow(geno_df_online) > 0) {
  hotspot_overlay_plot(
    gene_name = gene, gff_path = gff_path,
    annotations_df = ann_df_online, genotypes_df = geno_df_online
  )
}
} # }
```
