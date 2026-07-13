# Query genotypes for specific variant IDs (online)

Query genotypes for specific variant IDs (online)

## Usage

``` r
pg_query_genotypes(
  variant_ids,
  variant_id_col = "variant_id",
  variants_table = "variants",
  genotypes_table = "genotypes",
  meta_data = NULL
)
```

## Arguments

- variant_ids:

  A character vector of unique variant identifiers to look up (e.g.,
  `c("SNP_Chr05_75104557")`).

- variant_id_col:

  A character value specifying the primary key column name for variant
  IDs. Defaults to `"variant_id"`.

- variants_table:

  A character value specifying the name of the variants reference table.
  Defaults to `"variants"`.

- genotypes_table:

  A character value specifying the name of the genotype call data table.
  Defaults to `"genotypes"`.

- meta_data:

  A character vector of metadata columns to return alongside the
  genotype matrix. This can include any field from the \`variants\`
  table (e.g., \`"chrom"\`, \`"pos"\`, \`"ref"\`, \`"alt"\`,
  \`"variant_type"\`) as well as dynamically calculated metrics:
  \`"major_allele"\`, \`"minor_allele"\`, \`"major_allele_freq"\`, and
  \`"minor_allele_freq"\`. If \`NULL\` (default), all available database
  fields and all calculated metrics are returned.

## Value

A data frame in wide format (variants x samples) with the requested
metadata columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Define a list of target variant IDs
target_markers <- c("INDEL_Chr05_75104541", "SNP_Chr05_75104557")

# Query genotypes for these specific variant IDs from the online database,
# requesting specific metadata columns
pg_query_genotypes_result <- pg_query_genotypes(
  variant_ids = target_markers, 
  meta_data = c("chrom", "pos", "ref", "alt", "variant_type", 
                "minor_allele", "major_allele", "major_allele_freq", "minor_allele_freq")
)

# Print the head of the resulting data frame
print(head(pg_query_genotypes_result))
} # }
```
