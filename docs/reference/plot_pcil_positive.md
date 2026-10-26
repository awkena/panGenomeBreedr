# Plot PCIL Positive Introgressions

Visualizes the genomic regions of Pangenome-Characterized Introgression
Lines (PCIL) that were identified as positive for a given target. It
generates a plot for each target region, showing the introgression
blocks for all positive lines.

## Usage

``` r
plot_pcil_positive(pcil_positive_result)
```

## Arguments

- pcil_positive_result:

  A list object returned by \`fetch_pcil_positive\`, containing
  \`pcil_positive\` and \`regions\` data frames.

## Value

A list of ggplot objects, where each plot corresponds to a unique target
region and visualizes the introgression segments of the positive PCILs.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# 1. Connect to the package's bundled example database and fetch data
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)
con <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.
pcil_data <- fetch_pcil_data(con = con, connect_db_mode = "local")
selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")

variant_geno_sel <- fetch_genotypes_by_id(
  con = con,
  variant_ids = selection,
  connect_db_mode = "local"
)

# 2. Fetch relevant families
fam_results <- fetch_pcil_families_by_variant(
  con = con,
  selection = selection,
  pcil_data = pcil_data,
  connect_db_mode = "local"
)

# 3. Select positive lines
pcil_pos_pcv <- fetch_pcil_positive(
  pcil_data = pcil_data,
  variants_select_geno = variant_geno_sel,
  type = "position",
  sel = 15,
  available_ids = fam_results$pcil_summary[, c("sample_id", "selection")],
  result_pcil_families = fam_results,
  window = 0
)
#> Using +/- 0 bp window around positions

# 4. Generate and display the plot
positive_plots <- plot_pcil_positive(pcil_positive_result = pcil_pos_pcv)

if (length(positive_plots) > 0) {
  print(positive_plots[[1]])
}


disconnect_local_db(con)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
