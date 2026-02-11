# Identify SNP loci with potential genotype call errors.

Identify SNP loci with potential genotype call errors.

## Usage

``` r
geno_error(x, rp_row, dp_row, sep = ":", data_type = c("kasp", "agriplex"))
```

## Arguments

- x:

  A data frame of genotype calls where all columns are SNPs and samples
  as rows. Row names are unique sample names.

- rp_row:

  An integer or character value indication the row index or name of the
  recurrent or Parent 1.

- dp_row:

  An integer or character value indicating the row index or name of the
  donor or Parent 2.

- sep:

  A character used as separator for genotype calls, default is a colon.

- data_type:

  A character value indicating the data source; either \`kasp\` or
  \`agriplex\`.

## Value

A list object with the following components: 1) data frame of loci with
genotype call errors if present. 2) data frame of loci with genotype no
errors.

## Examples

``` r
# example code
# \donttest{
library(panGenomeBreedr)

# Marker data with errors in snp1 and snp4
dat <- data.frame(snp1 = c('C:A', 'A:A', 'C:A', 'C:T'),
                  snp2 = c('C:C', 'G:G', 'C:C', 'C:C'),
                  snp3 = c('C:T', 'C:C', 'C:T', 'C:T'),
                  snp4 = c('G:G', '-:-', 'G:-', 'G:A'),
                  snp5 = c('T:T', 'A:A', 'T:A', 'T:A'),
                  row.names = c('rp', 'dp', 'ind_1', 'ind_2'))

# Check for genotype call  error for each SNP
geno_mat <- geno_error(x = dat,
                  rp_row = 1,
                  dp_row = 2,
                  sep = ':',
                  data_type = 'kasp')
# }
```
