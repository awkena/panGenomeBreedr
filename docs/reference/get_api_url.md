# Get the current database API endpoint URL

Retrieves the active API endpoint. If no custom endpoint has been set by
the user, it safely defaults to the official Sorghum database API.

## Usage

``` r
get_api_url()
```

## Value

The character string of the active API URL.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the package
library(panGenomeBreedr)

# Get the current API endpoint URL
get_api_url()

} # }
```
