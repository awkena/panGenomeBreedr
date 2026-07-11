# Compute Linkage Disequilibrium (LD) Metrics (R2 and D')

Compute Linkage Disequilibrium (LD) Metrics (R2 and D')

## Usage

``` r
compute_LD(df, target_variant_ids = NULL, genotype_start_col = 11)
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
if (FALSE) { # \dontrun{
# Mode 1: Compute full pairwise matrix landscape

# Get genotype matrix
query_geno <- panGenomeBreedr::pg_query_db(
table_name = "genotypes",
chrom = "Chr03",
gene_name = "Sobic.003G421300",
start = 79037682,
end = 79039091
)

full_ld <- compute_LD(df = query_geno, target_variant_ids = NULL)
print(head(full_ld))

# Mode 2: Targeted calculation panel for KASP marker vetting
target_variants <- c("INDEL_Chr03_79037682", "INDEL_Chr03_79037750", "SNP_Chr03_79039022")
targeted_panel <- compute_LD(
  df = query_geno, 
  target_variant_ids = target_variants, 
  genotype_start_col = 11
)
} # }
```
