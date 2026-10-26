# Plot Best PCIL Candidate Lines

Creates a genome-wide overview plot for each target region, visualizing
the introgression patterns of the top-ranked Pangenome-Characterized
Introgression Lines (PCILs) identified by \`fetch_pcil_positive()\`.

## Usage

``` r
plot_pcil_best_lines(pcil_positive_result, pcil_data)
```

## Arguments

- pcil_positive_result:

  A list object returned by \`fetch_pcil_positive\`, which must contain
  the \`best_lines\` and \`regions\` data frames.

- pcil_data:

  A list object containing all required PCIL data tables, typically from
  \`fetch_pcil_data()\`. Must contain \`pcil_introgressions\`.

## Value

A list of ggplot objects. Each plot corresponds to a unique target
region and displays the introgression segments for the best candidate
lines across all chromosomes.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# 1. Connect to the package's bundled example database and find positive lines
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
fam_results <- fetch_pcil_families_by_variant(
  con = con,
  selection = selection,
  pcil_data = pcil_data,
  connect_db_mode = "local"
)

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

# 2. Plot the best candidate lines genome-wide
best_line_plots <- plot_pcil_best_lines(
  pcil_positive_result = pcil_pos_pcv,
  pcil_data = pcil_data
)

if (length(best_line_plots) > 0) {
  print(best_line_plots[[1]])
}


disconnect_local_db(con)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
