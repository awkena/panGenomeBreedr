# Color-code KASP genotype calls based on LGC genomics colors for HEX and FAM.

Color-code KASP genotype calls based on LGC genomics colors for HEX and
FAM.

## Usage

``` r
kasp_color(
  x,
  subset = "MasterPlate",
  sep = ":",
  geno_call = "Call",
  uncallable = "Uncallable",
  unused = "?",
  blank = "NTC",
  others = c("Missing", "Bad", "Dupe", "Over", "Short"),
  assign_cols = c(FAM = "blue", HEX = "gold", het = "forestgreen")
)
```

## Arguments

- x:

  A data frame of KASP genotype calls for one or multiple plates.

- subset:

  A character value indicating the column name for taking subsets of `x`
  for processing; default is the \`MasterPlate\`.

- sep:

  A character used as separator for genotype calls, default is a colon.

- geno_call:

  A character indicating the column name used for genotype calls in `x`.

- uncallable:

  A character indicating \`Uncallable\` genotype calls, if present.

- unused:

  A character indicating \`?\` genotype calls, if present.

- blank:

  A character value indicating \`No Template Controls (NTC)\` genotype
  calls.

- others:

  A character vector indicating other non-genotype calls in KASP
  genotype calls, if present. These may include \`'Missing', 'Bad',
  'Dupe'\`, \`'Over', 'Short'\`.

- assign_cols:

  A named character vector of \`length = 3\` for assigning colors to the
  FAM, HEX and heterozygous genotype groups.

## Value

A list object with subset unit as component data frames.

## Details

This is an experimental function. The default values of some of the
arguments in the function are based on LGC Genomics conventions
including the color codes for FAM and HEX fluorescence.

## Examples

``` r
# example code
library(panGenomeBreedr)
# \donttest{
dat1 <- kasp_color(x = panGenomeBreedr::kasp_dat,
                   subset = 'MasterPlate',
                   sep = ':',
                   geno_call = 'Call',
                   uncallable = 'Uncallable',
                   unused = '?',
                   blank = 'NTC')
#> Marker in Plate SE-24-0392_P01_d2 failed! Check genotype calls.
#> Marker in Plate SE-24-0392_P01_d1 failed! Check genotype calls.
#> Marker in Plate SE-24-0395_P01_d2 failed! Check genotype calls.
#> Marker in Plate SE-24-0395_P01_d1 failed! Check genotype calls.
#> Marker in Plate SE-24-0397_P01_d2 failed! Check genotype calls.
#> Marker in Plate SE-24-0397_P01_d1 failed! Check genotype calls.
# }
```
