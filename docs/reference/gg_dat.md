# Convert data matrix for genotypes to a long format data frame.

Convert data matrix for genotypes to a long format data frame.

## Usage

``` r
gg_dat(
  num_mat,
  map_file,
  map_pos = "pos",
  map_snp_ids = "snpid",
  map_chr = "chr"
)
```

## Arguments

- num_mat:

  A numeric or character data matrix with marker IDs as columns and
  sample IDs as row names.

- map_file:

  A data frame of map file consisting of SNP IDs and their chromosome
  numbers and positions as columns.

- map_pos:

  A character value indicating the column name for chromosome positions
  in map file.

- map_snp_ids:

  A character value indicating the column name for SNP IDs in
  `map_file`.

- map_chr:

  A character value indicating the column name for chromosome IDs in
  `map_file`.

## Value

A long format data frame object.

## Examples

``` r
# \donttest{
 # example code
library(panGenomeBreedr)
# Create a numeric matrix of genotype scores for 10 markers and 5 samples
num_dat <- matrix(c(rep(1, 10), rep(0, 10),
                    1, 1, 0.5, 1, 1, 1, 1, 1, 0, 1,
                    1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 0, 1, 1, 1, 1, 1, 1, 0.5 ),
                  byrow = TRUE, ncol = 10)

rownames(num_dat) <- c('rp', 'dp', paste0('bc1_', 1:3))
colnames(num_dat) <- paste0('S1', '_', c(floor(seq(1000, 10000, len = 8)),
                                         15000, 20000))

# Get map file by parsing SNP IDs
map_file <- parse_marker_ns(colnames(num_dat))

 # Convert num_geno to a long format data frame
 df <- gg_dat(num_mat = num_dat,
              map_file = map_file)
# }
```
