# Count the distribution of variant types in the database (online)

Count the distribution of variant types in the database (online)

## Usage

``` r
pg_count_variant_types(variants_table = "variants")
```

## Arguments

- variants_table:

  Character. The name of the table containing variant metadata. Defaults
  to "variants".

## Value

A data frame with two columns: `variant_type` and `n`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Count the distribution of variant types in the online database
pg_count_variant_types_result <- pg_count_variant_types()
print(pg_count_variant_types_result)
} # }

```
