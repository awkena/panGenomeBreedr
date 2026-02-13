# Get a summary of the number of samples per 96-well plate in a multi-plate KASP assay.

Get a summary of the number of samples per 96-well plate in a
multi-plate KASP assay.

## Usage

``` r
nsamples_plate(
  x,
  subset = "plates",
  snp_id = "SNPID",
  plate_id = "MasterPlate"
)
```

## Arguments

- x:

  A data frame of KASP genotype calls for one or multiple plates.

- subset:

  A character value indicating the column name for taking subsets of `x`
  for processing; default is the \`plates\`.

- snp_id:

  A character value indicating the column name for SNP IDs in `x`.

- plate_id:

  A character value indicating the column name for master plate having
  the same samples in `x`.

## Value

A list object with plates and a summary of number of samples per plate
as components.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)
dat1 <- panGenomeBreedr::beta_carotene
dat1 <- nsamples_plate(x = dat1,
                     subset = 'plates',
                     snp_id = 'SNPID',
                     plate_id = 'MasterPlate'
                    )$summ
# }
```
