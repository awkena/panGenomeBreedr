# Parse marker names of Hapmap format with a common pattern containing chromosome numbers and positions into a map file.

Parse marker names of Hapmap format with a common pattern containing
chromosome numbers and positions into a map file.

## Usage

``` r
parse_marker_ns(x, sep = "_", prefix = "S")
```

## Arguments

- x:

  A character vector containing the original marker names to be parsed.

- sep:

  A character value that serves as a unique separator between chromosome
  positions and other components of the marker name; default value is an
  underscore.

- prefix:

  A character value that represents a common pattern in marker names
  that precedes the chromosome number; the default value is \`S\`.

## Value

A data frame of map file consisting of the original marker names,
chromosome numbers and positions.

## Details

The marker names to be parsed into a map file must contain the
chromosome numbers and their positions with a common separator, as well
as a common pattern preceding the marker names. For instance, \`S1_101\`
and \`S2_102\` have \`S\` as the common pattern preceding the marker
names, \`1\` and \`2\` as the chromosome numbers, \`101\` and \`102\` as
the positions, and \`\_\` as the separator.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)

snps <- paste0('S', 1:10, '_', 101:110)
map_file <- parse_marker_ns(x = snps, sep = '_', prefix = 'S')
# }
```
