# Select polymorphic loci between two parents in a marker panel.

Select polymorphic loci between two parents in a marker panel.

## Usage

``` r
parent_poly(x, rp_row, dp_row, sep = ":", na_code = NA)
```

## Arguments

- x:

  A data matrix or frame with markers as columns and samples as rows.

- rp_row, dp_row:

  An integer or character value indicating the row index or name of
  Parent 1 and 2.

- sep:

  A character used as separator for genotype calls, default is a colon.

- na_code:

  A value indicating missing data in `x`.

## Value

A data matrix or frame object of polymorphic loci between two parental
lines.

## Examples

``` r
# example code

# Marker data
dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                  snp2 = c('C:C', 'C:C', 'C:C', 'C:C'),
                  snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
                  snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                  snp5 = c('T:T', 'T:T', 'T:A', 'T:A'),
                  row.names = c('rp', 'dp', 'art_het1', 'art_het2'))

# Find polymorphic loci
poly_loci <- parent_poly(x = dat,
                         rp_row = 1,
                         dp_row = 2,
                         sep = ':')
```
