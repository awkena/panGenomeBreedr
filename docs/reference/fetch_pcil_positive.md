# Fetch PCIL Positive Lines

Identifies and filters Pangenome-Characterized Introgression Lines
(PCIL) containing introgressions in a target gene, genomic region, or
position window. Supports multi-criteria ranking based on mean donor
fraction, haploblock size, and inbreeding coefficients (\$F\$).

## Usage

``` r
fetch_pcil_positive(
  pcil_data,
  variants_select_geno,
  type = c("gene", "region", "position"),
  sel = NULL,
  donor_thresh = 0.75,
  block_quantile = 0.75,
  F_quantile = 0.25,
  window = 0,
  available_ids = NULL,
  result_pcil_families = NULL,
  global_available_ids = NULL
)
```

## Arguments

- pcil_data:

  A list object containing all required PCIL data tables, typically from
  \`fetch_pcil_data()\`.

- variants_select_geno:

  The targets to search, in a form that depends on \`type\`:

  - `"gene"`: a character vector of gene IDs, or a data frame with a
    \`Region\` column of gene IDs.

  - `"region"`: a data frame with columns \`Region\` (a label), \`Chr\`,
    \`Start\` and \`End\`.

  - `"position"`: a data frame with \`variant_id\`, \`chrom\` and
    \`pos\` (as returned by \`fetch_genotypes_by_id()\`), or with
    \`Region\`, \`Chr\` and \`pos\`.

- type:

  Character string specifying the target type: \`"gene"\`, \`"region"\`
  or \`"position"\`.

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
  \`fetch_pcil_families_by_variant\`. When supplied, the search is
  restricted to sample IDs found in its \`pcil_summary\`.

- global_available_ids:

  Character vector of sample IDs. Restricts the search for every target
  to these lines only (for example, lines available in a seed bank).
  Applied to all targets, unlike \`available_ids\`, which is per target.
  If \`result_pcil_families\` is also supplied, only sample IDs present
  in both are used. Defaults to \`NULL\` (no restriction).

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

# 1. Connect to the package's bundled example database and fetch data/families
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)
con <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.
pcil_data <- fetch_pcil_data(con = con, connect_db_mode = "local")
selection <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")

results <- fetch_pcil_families_by_variant(
  con = con,
  selection = selection,
  pcil_data = pcil_data,
  connect_db_mode = "local"
)

# 2. Get genotypes for the selected variants
variant_geno_sel <- fetch_genotypes_by_id(
  con = con,
  variant_ids = selection,
  connect_db_mode = "local"
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

# 4. Select PCIL positives for an explicit genomic region
target_region <- data.frame(
  Region = "Chr03_target",
  Chr = "Chr03",
  Start = 79037800,
  End = 79037900
)
pcil_pos_region <- fetch_pcil_positive(
  pcil_data = pcil_data,
  variants_select_geno = target_region,
  type = "region"
)

disconnect_local_db(con)
#> Successfully disconnected from the local database. Memory cleared.
# }
```
