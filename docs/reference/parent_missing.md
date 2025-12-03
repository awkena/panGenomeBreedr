# Identify and subset loci with any parent missing genotype.

Identify and subset loci with any parent missing genotype.

## Usage

``` r
parent_missing(x, rp_row, dp_row, na_code = NA)
```

## Arguments

- x:

  A data matrix or frame with markers as columns and samples as rows.

- rp_row, dp_row:

  An integer or character value indicating the row index or name of
  Parent 1 and 2.

- na_code:

  A value indicating missing data in `x`.

## Value

A list object with the following components: 1) data frame of loci with
at least one parent genotype missing, if present. 2) data frame of loci
with all parent genotype present.

## Examples

``` r
# example code
library(panGenomeBreedr)

# Marker data
dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                  snp2 = c('C:C', NA, 'C:C', 'C:C'),
                  snp3 = c(NA, 'C:C', 'C:T', 'C:T'),
                  snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                  snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                  row.names = c('rp', 'dp', 'ind_1', 'ind_2'))

# Find loci with at least one missing parent genotype
par_miss <- parent_missing(x = dat,
                          rp_row = 1,
                          dp_row = 2)$par_missing
```
