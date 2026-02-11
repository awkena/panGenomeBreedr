# Get SNP or InDel alleles and possible genotypes from genotype calls in KASP data.

Get SNP or InDel alleles and possible genotypes from genotype calls in
KASP data.

## Usage

``` r
get_alleles(x, sep = ":", data_type = c("kasp", "agriplex"))
```

## Arguments

- x:

  A character vector of KASP genotype calls in one reaction plate.

- sep:

  A character used as separator for genotype calls, default is a colon.

- data_type:

  A character value indicating the data source; either \`kasp\` or
  \`agriplex\`.

## Value

A list object with \`length = 2\` consisting of marker alleles and
possible genotypes in `x`.

## Examples

``` r
# example code
# \donttest{
# Simulate a typical KASP genotype call
set.seed(123)
x <- sample(c("A:A", "A:-", "-:-", "Uncallable", "?"),
            size = 96,
            replace = TRUE)

# Assign NTC wells
x[c(88, 96)] <- 'NTC'

# Get alleles and expected genotypes
alleles <- get_alleles(x = x, data_type = 'kasp')
# }
```
