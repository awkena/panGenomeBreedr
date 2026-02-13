# Generate the prediction status of positive controls in a KASP assay, if present.

Generate the prediction status of positive controls in a KASP assay, if
present.

## Usage

``` r
pred_status(
  plate,
  geno_call = "Call",
  Group_id = "Group",
  blank = "NTC",
  Group_unknown = NULL
)
```

## Arguments

- plate:

  A data frame of KASP genotype calls for one plate.

- geno_call:

  A character indicating the column name used for genotype calls in
  `plate`.

- Group_id:

  A character value for the column ID indicating the predictions of the
  positive controls in `plate`.

- blank:

  A character value indicating \`No Template Controls (NTC)\` genotype
  calls.

- Group_unknown:

  A character value representing unknown genotype status for samples, if
  present. No genotype prediction can be made for such samples.

## Value

A data frame with the prediction status of each sample added as a
column.

## Details

The function recodes the prediction status of each sample as follows:
TRUE = prediction matches observed genotype call FALSE = prediction does
not match observed genotype call Unverified = Either observed genotype
call could not be made or expected genotype could not be made prior to
KASP genotyping or both. Blank = NTC wells

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)
# Get Plate 1

dat1 <- panGenomeBreedr::beta_carotene[1:96,]
dat1 <- pred_status(plate = dat1,
geno_call = 'Call',
Group_id = 'Group',
blank = 'NTC',
Group_unknown = '?')
# }
```
