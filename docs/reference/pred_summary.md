# Generate summary of prediction for positive controls in KASP genotype data, if present

Generate summary of prediction for positive controls in KASP genotype
data, if present

## Usage

``` r
pred_summary(
  x,
  snp_id = "SNPID",
  Group_id = NULL,
  blank = "NTC",
  Group_unknown = "?",
  geno_call = "Call",
  rate_out = FALSE
)
```

## Arguments

- x:

  A list object of KASP genotype calls processed by the \`kasp_color()\`
  function.

- snp_id:

  A character value indicating the column name for SNP IDs in `x`.

- Group_id:

  A character value for the column ID indicating the predictions of the
  positive controls in `x`.

- blank:

  A character value indicating \`No Template Controls (NTC)\` genotype
  calls.

- Group_unknown:

  A character value representing unverified expected genotype status for
  samples, if present. No genotype prediction can be made for such
  samples.

- geno_call:

  A character value indicating the column name for KASP genotype calls
  in `x`.

- rate_out:

  A logical value indicating whether to return raw counts or proportions
  of generated prediction summary.

## Value

A list object with plates and prediction summary as components.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)
dat1 <- panGenomeBreedr::beta_carotene
dat1 <- kasp_color(x = beta_carotene,
                   subset = 'plates',
                   sep = ':',
                   geno_call = 'Call',
                   uncallable = 'Uncallable',
                   unused = '?',
                   blank = 'NTC')

dat1 <- pred_summary(x = dat1,
                    snp_id = 'SNPID',
                    geno_call = 'Call',
                    Group_id = 'Group',
                    blank = 'NTC',
                    Group_unknown = '?')
dat1$summ
#>                          plate     snp_id false true unverified
#> 1 SE-24-1088_P01_d1_snpSB00800 snpSB00800     4    6         84
#> 2 SE-24-1088_P01_d2_snpSB00800 snpSB00800     2    6         86
#> 3 SE-24-1088_P01_d1_snpSB00803 snpSB00803     0   32         62
#> 4 SE-24-1088_P01_d2_snpSB00803 snpSB00803     0   32         62
#> 5 SE-24-1088_P01_d1_snpSB00804 snpSB00804     1   31         62
#> 6 SE-24-1088_P01_d2_snpSB00804 snpSB00804     1   31         62
#> 7 SE-24-1088_P01_d1_snpSB00805 snpSB00805    14   18         62
#> 8 SE-24-1088_P01_d2_snpSB00805 snpSB00805    14   18         62
# }
```
