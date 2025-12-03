# Calculate the proportion of recurrent parent background (RPP) fully recovered in backcross progenies.

Calculate the proportion of recurrent parent background (RPP) fully
recovered in backcross progenies.

## Usage

``` r
calc_rpp_bc(
  x,
  map_file,
  map_chr = "chr",
  map_pos = "pos",
  map_snp_ids = "snpid",
  rp = NULL,
  rp_num_code = 1,
  het_code = 0.5,
  na_code = NA,
  weighted = TRUE
)
```

## Arguments

- x:

  A numeric matrix of marker genotypes for backcross progenies and
  recurrent parent. Markers are columns and samples are rows.

- map_file:

  A data frame of map file consisting of marker IDs and their chromosome
  numbers and positions as columns.

- map_chr:

  A character value indicating the column name for chromosome IDs in
  `map_file`.

- map_pos:

  A character value indicating the column name for chromosome positions
  in map file.

- map_snp_ids:

  A character value indicating the column name for SNP IDs in
  `map_file`.

- rp:

  An integer or character value indicating the row index or sample ID
  for the recurrent parent. if \`NULL\` the sample on the first row of
  `x` is used as the recurrent parent.

- rp_num_code:

  A numeric value indicating the coding for the recurrent parent marker
  background in `x`.

- het_code:

  A numeric value indicating the coding for heterozygous marker
  background in `x`.

- na_code:

  A value indicating missing data in `x`.

- weighted:

  A logical value indicating whether RPP values should be weighted or
  not.

## Value

A data frame object comprising the RPP per chromosome and total RPP of
sample IDs.

## Details

The weighted RPP is computed as the sum of the product of relative
distance weighting of markers and genotype scores that match the
recurrent parent. Marker weights are obtained as half the normalized
ratio of the physical distance between marker \`i\` and marker \`i+1\`
to the total distance for all ordered markers on a chromosome.

The computation excludes heterozygous loci and only considers chromosome
regions fully recovered as recurrent parent genetic background.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)
# Create a numeric matrix of genotype scores for 10 markers and 5 samples
num_dat <- matrix(c(1, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 1, rep(0, 10),
                    1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
                    1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
                  byrow = TRUE, ncol = 10)

rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))

colnames(num_dat) <- c(paste0('S1', '_', c(floor(seq(1000, 10000, len = 5)))),
                            paste0('S2', '_', c(floor(seq(1000, 10000, len = 5)))))

# Get map file by parsing SNP IDs
map_file <- parse_marker_ns(colnames(num_dat))

# Calculate weighted RPP
rpp <- calc_rpp_bc(x = num_dat,
                      map_file = map_file,
                      map_chr = 'chr',
                      map_pos = 'pos',
                      map_snp_ids = 'snpid',
                      rp = 1,
                      rp_num_code = 1,
                      weighted = TRUE)
# }
```
