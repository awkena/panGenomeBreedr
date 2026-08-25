# Fetch Negative Control PCILs

Identifies and ranks negative control Pangenome-Characterized
Introgression Lines (PCILs) for a given set of positive lines. A
negative control is a line that is genetically similar (low
Identity-by-State distance) to a positive line but lacks the specific
introgression of interest.

## Usage

``` r
fetch_pcil_negative(
  pcil_data,
  pcil_positive_result,
  n_neg = 10,
  available_ids = NULL,
  result_pcil_families = NULL
)
```

## Arguments

- pcil_data:

  A list object containing PCIL data tables, typically from
  \`fetch_pcil_data()\`.

- pcil_positive_result:

  A list object returned by \`fetch_pcil_positive()\`, which must
  contain the \`best_lines\` and \`regions\` data frames.

- n_neg:

  Integer. The number of top-ranked negative control lines to return for
  each positive line. Defaults to \`10\`.

- available_ids:

  An optional data frame with \`sample_id\` and \`selection\` columns to
  restrict the search space for negative controls.

- result_pcil_families:

  An optional list object from \`fetch_pcil_families_by_variant()\` to
  globally restrict the search space.

## Value

A list containing:

- `pairs_best`: A data frame with the single best negative control for
  each positive line.

- `pairs_extended`: A data frame with up to \`n_neg\` ranked negative
  controls for each positive line.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# 1. Fetch data and establish positive lines
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

# 2. Fetch the top 10 negative control PCILs for the positive lines
pcil_neg_pcv <- fetch_pcil_negative(
  pcil_data = pcil_data,
  pcil_positive_result = pcil_pos_pcv,
  n_neg = 10
) 

print(head(pcil_neg_pcv$pairs_best))
#>                      SampleID_Positive  SampleID_Negative             Region
#> GMS_MN2025_132056   25ALM_BC1F3s1_2186  GMS_MN2025_132056 SNP_Chr03_79037855
#> GMS_MN2025_1320561  25ALM_BC1F3s1_0416  GMS_MN2025_132056 SNP_Chr03_79037855
#> 25ALM_BC1F3s1_0095  25ALM_BC1F3s1_1552 25ALM_BC1F3s1_0095 SNP_Chr03_79037855
#> GMS_MN2025_1320562   GMS_MN2025_114058  GMS_MN2025_132056 SNP_Chr03_79037855
#> 25ALM_BC1F3s1_00951 25ALM_BC1F3s1_0534 25ALM_BC1F3s1_0095 SNP_Chr03_79037855
#> GMS_MN2025_1320563  25ALM_BC1F3s1_1580  GMS_MN2025_132056 SNP_Chr03_79037855
#>                       Chr    Start      End   IBS_dis total_Mb_neg
#> GMS_MN2025_132056   Chr03 79037855 79037855 0.0870820        69.75
#> GMS_MN2025_1320561  Chr03 79037855 79037855 0.1251400        69.75
#> 25ALM_BC1F3s1_0095  Chr03 79037855 79037855 0.0910448        42.75
#> GMS_MN2025_1320562  Chr03 79037855 79037855 0.0737892        69.75
#> 25ALM_BC1F3s1_00951 Chr03 79037855 79037855 0.0917089        42.75
#> GMS_MN2025_1320563  Chr03 79037855 79037855 0.0906915        69.75
#>                     total_blocks_neg  F_neg Family Pair_Rank
#> GMS_MN2025_132056                  8 0.8332   SC49         1
#> GMS_MN2025_1320561                 8 0.8332   SC49         1
#> 25ALM_BC1F3s1_0095                 6 0.8706 SC1439         1
#> GMS_MN2025_1320562                 8 0.8332   SC49         1
#> 25ALM_BC1F3s1_00951                6 0.8706 SC1439         1
#> GMS_MN2025_1320563                 8 0.8332   SC49         1
# }
```
