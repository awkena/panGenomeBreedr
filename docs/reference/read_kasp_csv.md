# Read raw KASP results file (csv format) with one or multiple plates.

Read raw KASP results file (csv format) with one or multiple plates.

## Usage

``` r
read_kasp_csv(
  file,
  row_tags = c("Statistics", "DNA", "SNPs", "Scaling", "Data"),
  spacing = 2L,
  data_type = c("raw", "polished")
)
```

## Arguments

- file:

  A character value indicating the file name or path for the KASP
  results file (csv format).

- row_tags:

  A character vector for the ordered row tags for the components of the
  data in `file`.

- spacing:

  An integer value for specifying the number of empty rows between data
  components in `file`.

- data_type:

  A character value indicating the file type to import; currently
  supports either \`raw\` or \`polished\` KASP results file.

## Value

A list object of KASP results file for genotyping calls with FAM and HEX
coordinates.

## Examples

``` r
# example code
library(panGenomeBreedr)
# Read raw bulked KASP data into R
# \donttest{
path1 <-  system.file("extdata", "Genotyping_141.010_01.csv",
                      package = "panGenomeBreedr",
                      mustWork = TRUE)

file1 <- read_kasp_csv(file = path1, data_type = 'raw')
# Get KASP genotyping data for plotting
kasp_dat <- file1$Data
# }
```
