# Check whether progenies carry the donor segment at target loci

Tests the hypothesis that a trait locus was successfully introgressed.
For each progeny and each target locus, the markers in the locus are
classified as donor, heterozygous or recurrent parent calls. Each
progeny is then called at the locus and compared with the expected
genotype.

## Usage

``` r
introgression_check(
  x,
  map_file,
  snp_ids = "snpid",
  chr = "chr",
  chr_pos = "pos",
  parents,
  trait_pos,
  expected = c("donor", "het", "carrier"),
  threshold = 0.8,
  n_nearest = 3
)
```

## Arguments

- x:

  A numeric matrix of genotype scores with samples as rows and markers
  as columns, as returned by \`kasp_numeric()\`: \`1\` = recurrent
  parent allele, \`0\` = donor parent allele, \`0.5\` = heterozygous.
  Any other code (missing, error, monomorphic) is treated as not
  informative.

- map_file:

  A data frame of marker IDs, chromosome numbers and positions.

- snp_ids, chr, chr_pos:

  Column names in \`map_file\` for marker IDs, chromosome numbers and
  positions.

- parents:

  A character vector of length 2: the recurrent and donor parent row
  names in \`x\`. These rows are left out of the results.

- trait_pos:

  A named list of target loci. Each element is a named numeric vector
  with \`chr\` and either \`start\` and \`end\` (a region) or \`pos\` (a
  single position), as used by \`cross_qc_annotate()\`.

- expected:

  The genotype expected at the loci: \`"donor"\` (homozygous for the
  donor allele, e.g. NILs), \`"het"\` (heterozygous, e.g. BC1F1), or
  \`"carrier"\` (donor or heterozygous).

- threshold:

  Minimum proportion of informative markers that must agree for a call.

- n_nearest:

  Number of markers closest to the locus to use for a single position,
  or for a region that contains no markers.

## Value

A data frame with one row per progeny and locus: the locus, how the
markers were chosen (\`"region"\` or \`"nearest"\`), the number of
markers and informative markers, the proportions of donor, heterozygous
and recurrent calls, the call (\`Donor\`, \`Heterozygous\`,
\`Recurrent\`, \`Mixed\` or \`No informative markers\`) and whether it
matches the expected genotype.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)
num_dat <- matrix(c(rep(1, 6), rep(0, 6),
                    1, 1, 0, 0, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 0.5, 0.5, 1, 1),
                  byrow = TRUE, ncol = 6)
rownames(num_dat) <- c('rp', 'dp', paste0('bc_', 1:3))
colnames(num_dat) <- paste0('S1_', c(1000, 2000, 3000, 3500, 6000, 8000))
map_file <- parse_marker_ns(colnames(num_dat))

introgression_check(x = num_dat,
                    map_file = map_file,
                    parents = c('rp', 'dp'),
                    trait_pos = list(qtl1 = c(chr = 1, start = 2500, end = 4000)),
                    expected = 'donor')
#>   sample locus chr      region markers_from n_markers n_informative donor het
#> 1   bc_1  qtl1   1 2,500-4,000       region         2             2     1   0
#> 2   bc_2  qtl1   1 2,500-4,000       region         2             2     0   0
#> 3   bc_3  qtl1   1 2,500-4,000       region         2             2     0   1
#>   recurrent         call  pass
#> 1         0        Donor  TRUE
#> 2         1    Recurrent FALSE
#> 3         0 Heterozygous FALSE
# }
```
