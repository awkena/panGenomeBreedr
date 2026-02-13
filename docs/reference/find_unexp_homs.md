# Find loci with unexpected homozygous genotype calls for artificial heterozygotes.

Find loci with unexpected homozygous genotype calls for artificial
heterozygotes.

## Usage

``` r
find_unexp_homs(x, rp_row, dp_row, na_code = NA)
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
unexpected genotype calls artificial heterozygotes if present. 2) data
frame of loci with expected genotype calls for artificial heterozygotes.

## Details

Artificial heterozygotes are expected to show heterozygosity at
polymorphic loci between parents.Use this wrapper function to detect
loci which did not follow this prediction.

## Examples

``` r
# example code
library(panGenomeBreedr)
# Marker data
dat <- data.frame(snp1 = c('C:C', 'A:A', 'C:A', 'C:A'),
                  snp2 = c('C:C', '-:-', 'C:C', 'C:C'),
                  snp3 = c('T:T', 'C:C', 'C:T', 'C:T'),
                  snp4 = c('G:G', '-:-', 'G:-', 'G:G'),
                  snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                  row.names = c('rp', 'dp', 'art_het1', 'art_het2'))

# Check for unexpected homozygous genotypes
homs_unexp <- find_unexp_homs(x = dat,
                              rp_row = 1,
                              dp_row = 2)$geno_unexp
```
