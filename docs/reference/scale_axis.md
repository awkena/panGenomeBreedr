# Normalize FAM and HEX fluorescence values between 0 and 1

Normalize FAM and HEX fluorescence values between 0 and 1

## Usage

``` r
scale_axis(x)
```

## Arguments

- x:

  A numeric vector of FAM or HEX fluorescence values

## Value

A numeric vector of normalized values between 0 and 1

## Examples

``` r
# example code
library(panGenomeBreedr)
# Get Plate 1
dat1 <- panGenomeBreedr::kasp_dat[1:96,]
FAM_scaled <- scale_axis(dat1$X)
```
