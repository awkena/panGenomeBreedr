# Plot PCIL Positive and Negative Pairs

Creates a genome-wide overview plot for each positive-negative PCIL pair
identified by \`fetch_pcil_negative\`. Each plot visualizes the
introgression patterns for a positive line and its corresponding
negative controls.

## Usage

``` r
plot_pcil_pairs(pcil_neg_results, pcil_data)
```

## Arguments

- pcil_neg_results:

  A list object returned by \`fetch_pcil_negative\`, which must contain
  the \`pairs_extended\` data frame.

- pcil_data:

  A list object returned by \`fetch_pcil_data\`, containing all
  necessary PCIL data tables.

## Value

A list of ggplot objects. Each plot corresponds to a unique positive
PCIL and displays its introgression pattern alongside its ranked
negative control pairs.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# 1. Run the full pipeline to identify positive and negative pairs
pcil_data <- fetch_pcil_data(connect_db_mode = "online")
selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")

variant_geno_sel <- fetch_genotypes_by_id(variant_ids = selection, connect_db_mode = "online")
fam_results <- fetch_pcil_families_by_variant(
  selection = selection,
  pcil_data = pcil_data,
  connect_db_mode = "online"
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

pcil_neg_pcv <- fetch_pcil_negative(
  pcil_data = pcil_data,
  pcil_positive_result = pcil_pos_pcv,
  n_neg = 10
)

# 2. Generate side-by-side visual comparisons
pair_plots <- plot_pcil_pairs(
  pcil_neg_results = pcil_neg_pcv,
  pcil_data = pcil_data
)

if (length(pair_plots) > 0) {
  print(pair_plots[[1]])
}

# }
```
