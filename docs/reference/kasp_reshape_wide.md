# Reshape KASP data to wide format for same samples genotyped with multiple KASP markers.

Reshape KASP data to wide format for same samples genotyped with
multiple KASP markers.

## Usage

``` r
kasp_reshape_wide(
  x,
  subset = "MasterPlate",
  snp_id = "SNPID",
  geno_call = "Call",
  idvar = "SubjectID",
  blank = "NTC"
)
```

## Arguments

- x:

  A data frame of KASP genotype calls for one or multiple plates.

- subset:

  A character value indicating the column name for master plate having
  the same samples in `x`.

- snp_id:

  A character value indicating the column name for SNP IDs in `x`.

- geno_call:

  A character value indicating the column name for KASP genotype calls
  in `x`.

- idvar:

  A character value indicating the column name for unique subject IDs of
  samples in `x`.

- blank:

  A character value indicating \`No Template Controls (NTC)\` genotype
  calls.

## Value

A list object with reshaped data for master plates that have the same
subject IDs as the components. The components are data frames with
Column 1 as the subject IDs and the rest of the columns as the number of
KASP markers assayed for each master plate.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)
dat1 <- panGenomeBreedr::beta_carotene
plate_wide <- kasp_reshape_wide(x = dat1,
                                subset = 'MasterPlate',
                                snp_id = 'SNPID',
                                geno_call = 'Call',
                                idvar = "SubjectID",
                                blank = 'NTC')
# }
```
