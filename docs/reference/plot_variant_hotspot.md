# Plot genomic variant hotspots (SNPs and INDELs) extracted from the database (local)

Visualizes variant distributions across a genomic region. Variants are
mapped by position, with impact levels determining their vertical
placement (lollipop style) and variant types mapped to shapes.

## Usage

``` r
plot_variant_hotspot(
  var_df,
  chrom_col = "chrom.x",
  pos_col = "pos.x",
  start_pos = NULL,
  end_pos = NULL
)
```

## Arguments

- var_df:

  A data frame containing variant data. Required columns: `variant_type`
  (character: "SNP", "INDEL"), `variant_id`, and `impact` (character:
  "MODIFIER", "LOW", "MODERATE", "HIGH").

- chrom_col:

  A character string specifying the chromosome column.

- pos_col:

  A character string specifying the variant position column.

- start_pos:

  A numeric integer for the region start. Default is the minimum
  position in var_df.

- end_pos:

  A numeric integer for the region end. Default is the maximum position
  in var_df.

## Value

A `ggplot` object representing the variant hotspots.

## Examples

``` r
# \donttest{
library(panGenomeBreedr)
# Extract variants within the candidate gene: Sobic.005G213600
pg_gt_region <- pg_query_db(table_name = "variants",
                            chrom = "Chr05",
                            start = 75104537,
                            end = 75106403
                            )

#' # Extract annotations for the candidate gene
pg_annota_region <- pg_query_db(table_name = "annotations",
                                chrom = "Chr05",
                                start = 75104537,
                                end = 75106403,
                                gene_name = "Sobic.005G213600"
                                )

var_df <- merge(pg_annota_region, pg_gt_region, by = "variant_id")

plot_variant_hotspot(var_df)

# }
```
