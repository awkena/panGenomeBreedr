# Identify lines that possess favorable alleles for target loci using trait predictive markers.

Identify lines that possess favorable alleles for target loci using
trait predictive markers.

## Usage

``` r
foreground_select(
  geno_data,
  fore_marker_info,
  fore_marker_col,
  fav_allele_col,
  alt_allele_col,
  select_type = c("homo", "hetero", "both"),
  sep = ":"
)
```

## Arguments

- geno_data:

  A data frame or matrix of marker genotype data with lines as rows and
  markers as columns.

- fore_marker_info:

  A data frame of meta data for the trait-predictive markers consisting
  of marker names, favorable and alternate alleles as variables.

- fore_marker_col:

  A character value specifying the column name for trait-predictive
  markers in `fore_marker_info`.

- fav_allele_col:

  A character value specifying the column name for favorable allele in
  `fore_marker_info`.

- alt_allele_col:

  A character value specifying the column name for alternate allele in
  `fore_marker_info`.

- select_type:

  A character value of one of three options: \`homo\` to select lines
  that are homozygous for the favorable allele; \`hetero\` to select
  lines that are heterozygous; and \`both\` to select both favorable
  homozygotes and heterozygotes.

- sep:

  A character value representing the separator used for genotype
  calling.

## Value

A binary data frame conversion of the original marker genotype data for
presence(1s) and absence(0s) of favorable alleles of target loci.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)

# Marker genotype data
geno <- data.frame(SNP1 = c("A:A", "A:G", "G:G", "A:A"),
                   SNP2 = c("C:C", "C:T", "T:T", "C:T"),
                   SNP3 = c("G:G", "G:G", "A:G", "A:A"),
                   row.names = c("Line1", "Line2", "Line3", "Line4"))

# Trait predictive markers meta data
marker_info <- data.frame(qtl_markers = paste0('SNP', 1:3),
                          fav_alleles = c('A', 'C', 'G'),
                          alt_alleles = c('G', 'T', 'A'))

# Select lines where genotype is either homozygous for favorable alleles at all loci
foreground_select(geno_data = geno,
                  fore_marker_info = marker_info,
                  fore_marker_col = 'qtl_markers',
                  fav_allele_col = 'fav_alleles',
                  alt_allele_col = 'alt_alleles',
                  select_type = "homo")
#>       SNP1 SNP2 SNP3
#> Line1    1    1    1
#> Line2    0    0    1
#> Line3    0    0    0
#> Line4    1    0    0

# }
```
