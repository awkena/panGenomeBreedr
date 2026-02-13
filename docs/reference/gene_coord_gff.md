# Get the genomic range of a candidate gene using the Sobic ID from a gff file.

Get the genomic range of a candidate gene using the Sobic ID from a gff
file.

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

A list object of three components indicating the chromosome, start and
end coordinates of candidate gene.

## Examples

``` r
# \donttest{
# example code

# Path to GFF3 file
gff_path <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
gene_coord_gff(gene_name = "Sobic.005G213600",
               gff_path = gff_path)
#> Warning: downloaded length 0 != reported length 454
#> Warning: cannot open URL 'https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz': HTTP status was '503 Service Unavailable'
#> Error in utils::download.file(gff_path, destfile = file.path(tempdir(),     basename(gff_path))): cannot open URL 'https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz'

# }
```
