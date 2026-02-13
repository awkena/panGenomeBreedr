# Order marker IDs based on their chromosome numbers and positions in ascending order.

Order marker IDs based on their chromosome numbers and positions in
ascending order.

## Usage

``` r
order_markers(x, chr_col = "chr", pos_col = "pos")
```

## Arguments

- x:

  A data frame object for a map file or marker data file.

- chr_col:

  A character value indicating the column name for chromosome IDs in
  `x`.

- pos_col:

  A character value indicating the column name for chromosome positions
  in `x`.

## Value

A data frame object of the same dimension as `x`.

## Examples

``` r
# example code
map_file <- data.frame(snpid = paste0('S', rep(1:2, 5), '_', 1001:1005),
                       chr = rep(1:2, 5),
                       pos = rep(1001:1005, 2))

# Order map file
map_file <- order_markers(x = map_file,
                          chr_col = 'chr',
                          pos_col = 'pos')
```
