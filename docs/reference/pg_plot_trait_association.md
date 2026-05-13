# Generic Trait Association Audit (Cloud/API Version)

This function performs an association analysis between a specific
genomic variant and a phenotypic trait. It pulls genotypes from the AWS
API and merges them with your local phenotypic data for statistical
testing and plotting.

## Usage

``` r
pg_plot_trait_association(
  variant_id,
  pheno_df,
  trait_col,
  trait_type = c("qualitative", "quantitative")
)
```

## Arguments

- variant_id:

  Character. The specific SNP or INDEL ID to analyze.

- pheno_df:

  Data frame containing local phenotypic data. The first column must
  contain accession identifiers (PI numbers).

- trait_col:

  Character. The name of the column in `pheno_df` containing the trait
  scores.

- trait_type:

  Character. Either "qualitative" or "quantitative". Determines the
  statistical test and plot type.

## Value

A ggplot object displaying the association and statistical summary.
