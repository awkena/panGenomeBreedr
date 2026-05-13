# Get the genomic range of a candidate gene using the Sobic ID from a GFF file

This function parses a GFF3 file to extract the chromosome, start, and
end coordinates for a specific gene ID. It supports local files,
GZ-compressed files, and direct URLs.

## Usage

``` r
pg_gene_coords(gene_name, gff_path)
```

## Arguments

- gene_name:

  A character value indicating the Sobic ID of the candidate gene.

- gff_path:

  A character value specifying the path to the GFF3 file. URL paths and
  GZ-compressed files are supported.

## Value

A list containing the chromosome, start, and end coordinates.

## Examples

``` r
if (FALSE) { # \dontrun{
library(panGenomeBreedr)

# Define path to a sorghum GFF3 file (v5.1)
gff_url <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"

# Retrieve coordinates for a candidate gene
coords <- pg_gene_coords(gene_name = "Sobic.005G213600", gff_path = gff_url)
print(coords)
} # }

```
