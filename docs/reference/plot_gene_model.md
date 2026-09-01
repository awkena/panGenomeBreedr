# Plot the gene model using genomic coordinates (local)

Visualizes the gene structure with adjusted component heights, a cleaned
legend, and a strand orientation arrow. Introns are implicitly
represented by the exposed transcript backbone. Each transcript is drawn
on its own row below the gene backbone (labeled on the y-axis), so
alternative transcripts for the same gene are shown as separate tracks
instead of being collapsed on top of one another.

## Usage

``` r
plot_gene_model(gene_df, row_height = 0.6)
```

## Arguments

- gene_df:

  A data frame containing the genomic coordinates, expected to be the
  output from `gene_coord_gff`.

- row_height:

  A numeric value controlling the vertical spacing between transcript
  rows. Default is 0.6.

## Value

A `ggplot` object representing the structural gene model.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)
# Path to GFF3 file
gff_path <- "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
gene_features <- gene_coord_gff(gene_name = "Sobic.005G213600",
                                gff_path = gff_path)
# Plot gene model for Sobic.005G213600
plot_gene_model(gene_df = gene_features)

# }
```
