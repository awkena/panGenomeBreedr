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
# Mode 1: Compute full pairwise matrix landscape

# Get genotype matrix from the online database
query_geno <- panGenomeBreedr::fetch_table_region(
table_name = "genotypes",
chrom = "Chr03",
start = 79037682,
end = 79039091,
connect_db_mode = 'online'
)

full_ld <- calculate_LD(df = query_geno, target_variant_ids = NULL, genotype_start_col = 11)
print(head(full_ld))
#>            variant_1 position_1 variant_type_1          variant_2 position_2
#> 1 SNP_Chr03_79037693   79037693            SNP SNP_Chr03_79037699   79037699
#> 2 SNP_Chr03_79037693   79037693            SNP SNP_Chr03_79037706   79037706
#> 3 SNP_Chr03_79037693   79037693            SNP SNP_Chr03_79037716   79037716
#> 4 SNP_Chr03_79037693   79037693            SNP SNP_Chr03_79037855   79037855
#> 5 SNP_Chr03_79037693   79037693            SNP SNP_Chr03_79037876   79037876
#> 6 SNP_Chr03_79037693   79037693            SNP SNP_Chr03_79037944   79037944
#>   variant_type_2 distance_bp    R2 D_prime
#> 1            SNP           6 2e-05       1
#> 2            SNP          13 0e+00       1
#> 3            SNP          23 0e+00       1
#> 4            SNP         162 5e-05       1
#> 5            SNP         183 1e-05       1
#> 6            SNP         251 5e-05       1

# Mode 2: Targeted calculation panel for KASP marker vetting
target_variants <- c("INDEL_Chr03_79037682", "INDEL_Chr03_79037750", "SNP_Chr03_79039022")
targeted_panel <- calculate_LD(
  df = query_geno, 
  target_variant_ids = target_variants, 
  genotype_start_col = 11
)
print(head(targeted_panel))
#>              variant_1 position_1 variant_type_1          variant_2 position_2
#> 1 INDEL_Chr03_79037682   79037682          INDEL SNP_Chr03_79037693   79037693
#> 2 INDEL_Chr03_79037682   79037682          INDEL SNP_Chr03_79037699   79037699
#> 3 INDEL_Chr03_79037682   79037682          INDEL SNP_Chr03_79037706   79037706
#> 4 INDEL_Chr03_79037682   79037682          INDEL SNP_Chr03_79037716   79037716
#> 5 INDEL_Chr03_79037682   79037682          INDEL SNP_Chr03_79037855   79037855
#> 6 INDEL_Chr03_79037682   79037682          INDEL SNP_Chr03_79037876   79037876
#>   variant_type_2 distance_bp      R2 D_prime
#> 1            SNP          11 0.00001 1.00000
#> 2            SNP          17 0.00048 1.00000
#> 3            SNP          24 0.00001 1.00000
#> 4            SNP          34 0.00002 1.00000
#> 5            SNP         173 0.00010 0.02111
#> 6            SNP         194 0.00035 1.00000
# }
```
