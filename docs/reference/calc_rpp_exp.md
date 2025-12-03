# Compute theoretical RPP values for any specified backcross generation.

Compute theoretical RPP values for any specified backcross generation.

## Usage

``` r
calc_rpp_exp(bc_gen = 1, rpp2n = FALSE)
```

## Arguments

- bc_gen:

  A positive integer value indicating the backcross generation.

- rpp2n:

  A logical value indicating whether to compute theoretical RPP values
  for the F1 all BC generations up to `bc_gen`.

## Value

A named vector of theoretical RPP values.

## Examples

``` r
# example code
# Calculate theoretical RPP values up to BC5
rpp_values <- calc_rpp_exp(bc_gen = 5, rpp2n = TRUE)
```
