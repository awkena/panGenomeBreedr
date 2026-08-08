# Get the genomic range of a candidate gene and its features using the Sobic ID from a GFF file (local)

Get the genomic range of a candidate gene and its features using the
Sobic ID from a GFF file (local)

## Usage

``` r
gene_coord_gff(gene_name, gff_path)
```

## Arguments

- gene_name:

  A character value indicating the Sobic ID of candidate gene.

- gff_path:

  A character value specifying the path to gff3 file. GZ compressed
  files and URL file paths are acceptable.

## Value

A data frame containing the genomic coordinates (chromosome, start, end,
strand) of the gene, transcripts, 5'-UTR, CDS, and 3'-UTR.

## Examples

``` r
# \donttest{
# example code
library(panGenomeBreedr)
# Path to GFF3 file
gff_path <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
gene_features <- gene_coord_gff(gene_name = "Sobic.005G213600",
                                gff_path = gff_path)
head(gene_features)
#>                        ID         Feature Chromosome    Start      End Strand
#> 1 Sobic.005G213600.1.v5.1            mRNA      Chr05 75104537 75106403      -
#> 2 Sobic.005G213600.1.v5.1 three_prime_UTR      Chr05 75104537 75104865      -
#> 3 Sobic.005G213600.1.v5.1             CDS      Chr05 75104866 75106168      -
#> 4 Sobic.005G213600.1.v5.1             CDS      Chr05 75106279 75106334      -
#> 5 Sobic.005G213600.1.v5.1  five_prime_UTR      Chr05 75106335 75106403      -
#> 6   Sobic.005G213600.v5.1            gene      Chr05 75104537 75106403      -
# }
```
