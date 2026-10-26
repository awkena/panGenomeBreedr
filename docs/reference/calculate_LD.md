# Compute Linkage Disequilibrium (LD) Metrics (R2 and D')

Compute Linkage Disequilibrium (LD) Metrics (R2 and D')

## Usage

``` r
calculate_LD(df, target_variant_ids = NULL, genotype_start_col = 11)
```

## Arguments

- df:

  A data frame containing the genotype matrix. Must include a
  \`variant_id\` column and a genomic position column (e.g., 'pos',
  'position').

- target_variant_ids:

  A character vector of target variant IDs. If provided, the function
  runs in "Targeted Mode". If \`NULL\` (default), it runs in "All-vs-All
  Mode".

- genotype_start_col:

  An integer specifying the column index where phased genotype calls
  (e.g., "0\|0", "1\|0") begin. Defaults to 7.

## Value

A long-format data frame with the following columns:

- variant_1:

  ID of the first variant in the pair.

- variant1_position:

  Genomic position of the first variant.

- variant1_type:

  Type of the first variant (SNP or INDEL).

- variant_2:

  ID of the second variant in the pair.

- variant2_position:

  Genomic position of the second variant.

- variant2_type:

  Type of the second variant (SNP or INDEL).

- distance_bp:

  Absolute physical distance in base pairs between the variants.

- R2:

  The squared correlation coefficient (\$R^2\$).

- D_prime:

  The normalized linkage disequilibrium coefficient (\$D'\$).

## Examples

``` r
# \donttest{
library(panGenomeBreedr)

# Locate the package example database folder
my_db_folder <- system.file("extdata", "pangenome_scale_db",
                           package = "panGenomeBreedr",
                           mustWork = TRUE)
con_local <- connect_local_db(folder_path = my_db_folder)
#> Successfully connected to the local offline database! Pangenome-scale database  mounted safely.

# Get genotype matrix for the region from the local database
query_geno <- fetch_table_region(
  con = con_local,
  table_name = "genotypes",
  chrom = "Chr03",
  start = 79037682,
  end = 79039091
)

# Mode 1: Compute full pairwise matrix landscape
full_ld <- calculate_LD(df = query_geno, target_variant_ids = NULL, genotype_start_col = 11)
print(head(full_ld))
#>              variant_1 position_1 variant_type_1            variant_2
#> 1 INDEL_Chr03_79037682   79037682          INDEL   SNP_Chr03_79037693
#> 2 INDEL_Chr03_79037682   79037682          INDEL   SNP_Chr03_79037699
#> 3 INDEL_Chr03_79037682   79037682          INDEL   SNP_Chr03_79037706
#> 4 INDEL_Chr03_79037682   79037682          INDEL   SNP_Chr03_79037716
#> 5 INDEL_Chr03_79037682   79037682          INDEL INDEL_Chr03_79037750
#> 6 INDEL_Chr03_79037682   79037682          INDEL   SNP_Chr03_79037855
#>   position_2 variant_type_2 distance_bp      R2 D_prime
#> 1   79037693            SNP          11 0.00001 1.00000
#> 2   79037699            SNP          17 0.00048 1.00000
#> 3   79037706            SNP          24 0.00001 1.00000
#> 4   79037716            SNP          34 0.00002 1.00000
#> 5   79037750          INDEL          68 0.00010 1.00000
#> 6   79037855            SNP         173 0.00010 0.02111

# Mode 2: Targeted calculation panel for KASP marker vetting
target_variants <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")
targeted_panel <- calculate_LD(
  df = query_geno,
  target_variant_ids = target_variants,
  genotype_start_col = 11
)
print(head(targeted_panel))
#>              variant_1 position_1 variant_type_1            variant_2
#> 1 INDEL_Chr03_79037889   79037889          INDEL INDEL_Chr03_79037682
#> 2 INDEL_Chr03_79037889   79037889          INDEL   SNP_Chr03_79037693
#> 3 INDEL_Chr03_79037889   79037889          INDEL   SNP_Chr03_79037699
#> 4 INDEL_Chr03_79037889   79037889          INDEL   SNP_Chr03_79037706
#> 5 INDEL_Chr03_79037889   79037889          INDEL   SNP_Chr03_79037716
#> 6 INDEL_Chr03_79037889   79037889          INDEL INDEL_Chr03_79037750
#>   position_2 variant_type_2 distance_bp      R2 D_prime
#> 1   79037682          INDEL         207 0.00007 0.02088
#> 2   79037693            SNP         196 0.00006 1.00000
#> 3   79037699            SNP         190 0.00283 1.00000
#> 4   79037706            SNP         183 0.00006 1.00000
#> 5   79037716            SNP         173 0.00012 1.00000
#> 6   79037750          INDEL         139 0.00002 0.01667

disconnect_local_db(con_local)
#> Successfully disconnected from the local database. Memory cleared.
# }

if (FALSE) { # \dontrun{
# To query the public online resource instead:
query_geno_online <- fetch_table_region(
  table_name = "genotypes",
  chrom = "Chr03",
  start = 79037682,
  end = 79039091,
  connect_db_mode = 'online'
)
full_ld_online <- calculate_LD(df = query_geno_online, target_variant_ids = NULL,
                                genotype_start_col = 11)
print(head(full_ld_online))
} # }
```
