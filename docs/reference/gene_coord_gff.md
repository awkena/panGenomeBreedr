# Get the genomic range of a candidate gene from a GFF file

This function parses a GFF3 file (local or remote) to find the genomic
coordinates (chromosome, start, end) for a specified gene ID.

## Usage

``` r
gene_coord_gff(gene_name, gff_path)
```

## Arguments

- gene_name:

  A character value indicating the Sobic ID of the candidate gene.

- gff_path:

  A character value specifying the path to the GFF3 file. URL paths and
  GZ-compressed files are supported.

## Value

A list object of three components indicating the chromosome, start and
end coordinates of candidate gene.

## Examples

``` r
# \donttest{
# example code

# Define path to a remote sorghum GFF3 file (v5.1)
gff_path <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"

# Retrieve coordinates for a candidate gene
gene_coord_gff(gene_name = "Sobic.005G213600",
               gff_path = gff_path)
#> $chrom
#> [1] "Chr05"
#> 
#> $start
#> [1] 75104537
#> 
#> $end
#> [1] 75106403
#> 

# }
```
