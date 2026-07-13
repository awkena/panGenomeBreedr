# Filter genotypes by multiple sample metadata criteria

Filter genotypes by multiple sample metadata criteria

## Usage

``` r
pg_query_geno_by_meta(genotype_matrix, genotype_start_col, filters)
```

## Arguments

- genotype_matrix:

  A data frame or matrix where rows are variants and columns include
  variant metadata followed by sample genotypes.

- genotype_start_col:

  An integer specifying the column index where sample genotype calls
  begins.

- filters:

  A named list of metadata criteria to filter samples by. Names must
  match metadata columns, and values are the allowed entries. Example:
  `list(population = "Gates", countryorigin = c("Ethiopia", "Ghana", "Togo"))`

## Value

A data frame containing the filtered genotype matrix with recalculated
allele metrics for the sub-population. Returns an empty data frame if no
samples match all combined criteria.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Define filtering criteria
my_filters <- list(
  population = "Gates",
  countryorigin = c("Ethiopia", "Ghana", "Togo")
)

# Extract the full genotype matrix first for a genomic region
genotype_data_region <- pg_query_db(
  table_name = "genotypes",
  chrom = "Chr05",
  start = 75104537,
  end = 75106403
)

# Define the genotype_start_col.
genotype_start_col_val <- 11

# Get filtered genotypes matrix
filtered_genotypes <- pg_query_geno_by_meta(
  genotype_matrix = genotype_data_region,
  genotype_start_col = genotype_start_col_val,
  filters = my_filters
)
print(head(filtered_genotypes))
} # }
```
