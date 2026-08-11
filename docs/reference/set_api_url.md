# Set the API endpoint URL for your private database

Set the API endpoint URL for your private database

## Usage

``` r
set_api_url(url)
```

## Arguments

- url:

  Character. The full URL to your hosted API endpoint (e.g.,
  "http://rice-genomics.myuniversity.edu:8000").

## Examples

``` r
# \donttest{
# Load the package
library(panGenomeBreedr)

# Set a custom API endpoint URL
set_api_url("http://132.145.61.28:8000")
#> panGenomeBreedr API endpoint successfully set to: http://132.145.61.28:8000

# }
```
