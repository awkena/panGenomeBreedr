# Format marker names to comply with the Hapmap convention.

Format marker names to comply with the Hapmap convention.

## Usage

``` r
hapmap_ns_fmt(
  x,
  map_file,
  snpid_col = "SNP_name",
  chr_col = "chr",
  pos_col = "pos",
  kasp_prefix = "snpSB",
  scaffold_prefix = "Sbv"
)
```

## Arguments

- x:

  A character vector of marker IDs to be formated.

- map_file:

  A data frame of map file consisting of marker IDs and their chromosome
  numbers and positions as columns.

- snpid_col:

  A character value indicating the column name for marker IDs in
  `map_file`.

- chr_col:

  A character value indicating the column name for chromosome IDs in
  `map_file`.

- pos_col:

  A character value indicating the column name for chromosome positions
  in `map_file`.

- kasp_prefix:

  A character value indicating the KASP marker prefix for the crop
  species.

- scaffold_prefix:

  A character value indicating the KASP marker prefix for the crop
  species.

## Value

A character vector of markers names formated to comply with the Hapmap
format.

## Details

This is an experimental function, so use it with caution. It can be used
to create marker IDs that can easily be parsed into a map file using the
\`parse_marker_ns()\` function. The function may not format all inputted
marker names in `x`.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)

# Marker IDs
snps <- c('snpSB00072', 'snpSB00106', 'snpSB00109', 'Sbv3.1_01_68028666I',
          'Sbv3.1_02_67884158W')

# Map file for SNPs
map_file <- data.frame(snpid = snps,
                       chr = c(2, 5, 5, 1, 2),
                       pos = c(61811307, 838874, 1730282, NA, NA))

# Format marker IDs to hapmap format
ns_new <- hapmap_ns_fmt(x = snps,
                        map_file = map_file,
                        snpid_col = 'snpid',
                        chr_col = 'chr',
                        pos_col = 'pos',
                        kasp_prefix = 'snpSB',
                        scaffold_prefix = 'Sbv')
# }
```
