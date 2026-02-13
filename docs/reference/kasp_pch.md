# Generate pch characters for cluster plots of KASP genotype calls.

Generate pch characters for cluster plots of KASP genotype calls.

## Usage

``` r
kasp_pch(
  x,
  sep = ":",
  blank = "NTC",
  uncallable = "Uncallable",
  unused = "?",
  others = c("Missing", "Bad", "Dupe", "Over", "Short")
)
```

## Arguments

- x:

  A character vector of KASP genotype calls in one reaction plate.

- sep:

  A character used as separator for genotype calls, default is a colon.

- blank:

  A character value indicating \`No Template Controls (NTC)\` genotype
  calls.

- uncallable:

  A character indicating \`Uncallable\` genotype calls, if present.

- unused:

  A character indicating \`?\` genotype calls, if present.

- others:

  A character vector indicating other non-genotype calls in KASP
  genotype calls, if present. These may include \`'Missing', 'Bad',
  'Dupe'\`, \`'Over', 'Short'\`.

## Value

A \`96 x 1\` data frame of pch values for possible genotypes in each
KASP reaction plate.

## Examples

``` r
# example code
# \donttest{
x <- panGenomeBreedr::kasp_dat$Call[1:96]
geno_pch <- kasp_pch(x = x)
# }
```
