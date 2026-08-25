# Fetch PCIL Positive Lines

Identifies and filters Pangenome-Characterized Introgression Lines
(PCIL) containing introgressions in a target gene or position window.
Supports multi-criteria ranking based on mean donor fraction, haploblock
size, and inbreeding coefficients (\$F\$).

## Usage

``` r
fetch_pcil_positive(
  pcil_data,
  variants_select_geno,
  type = c("gene", "position"),
  sel = NULL,
  donor_thresh = 0.75,
  block_quantile = 0.75,
  F_quantile = 0.25,
  window = 0,
  available_ids = NULL,
  result_pcil_families = NULL
)
```

## Arguments

- pcil_data:

  A list object containing all required PCIL data tables, typically from
  \`fetch_pcil_data()\`.

- variants_select_geno:

  A character vector (gene IDs) or a data frame containing position
  definitions.

- type:

  Character string specifying the target type: \`"gene"\` or
  \`"position"\`.

- sel:

  Integer. Number of top-ranked candidate lines to select per region.
  Defaults to \`NULL\` (returns all).

- donor_thresh:

  Numeric. Minimum mean donor fraction required. Defaults to \`0.75\`.

- block_quantile:

  Numeric. Upper quantile threshold to exclude large introgressions.
  Defaults to \`0.75\`.

- F_quantile:

  Numeric. Lower quantile threshold for the inbreeding coefficient
  (\$F\$). Defaults to \`0.25\`.

- window:

  Integer. Symmetric base pair window (\$+/-\$) around targets when
  \`type = "position"\`. Defaults to \`100\`.

- available_ids:

  Data frame with columns \`sample_id\` and \`selection\` to restrict
  line searches per target.

- result_pcil_families:

  list object of results obtained from
  \`fetch_pcil_families_by_variant\`.

## Value

A list containing:

- `pcil_positive`: Data frame of all matching positive introgression
  lines merged with genome-wide stats.

- `best_lines`: (Optional) Top-ranked candidate lines per region if
  `sel` is specified.

- `regions`: Data frame of normalized genomic target coordinates
  evaluated.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# 1. Fetch data and families
pcil_data <- fetch_pcil_data(connect_db_mode = 'online')
selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")

results <- fetch_pcil_families_by_variant(
  selection = selection,
  pcil_data = pcil_data,
  connect_db_mode = 'online'
)

# 2. Get genotypes for the selected variants
variant_geno_sel <- fetch_genotypes_by_id(
  variant_ids = selection, 
  connect_db_mode = 'online'
)

# 3. Select PCIL positives by variant positions
pcil_pos_pcv <- fetch_pcil_positive(
  pcil_data = pcil_data,
  variants_select_geno = variant_geno_sel,
  type = "position",
  sel = 15,
  available_ids = results$pcil_summary[c("sample_id", "selection")],
  result_pcil_families = results,
  window = 0
)
#> Using +/- 0 bp window around positions
# }
```
