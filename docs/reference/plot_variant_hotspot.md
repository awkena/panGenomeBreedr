# Plot genomic variant hotspots (SNPs and INDELs) extracted from the database (local)

Visualizes variant distributions across a genomic region. Variants are
mapped by position, with impact levels determining their vertical
placement (lollipop style), variant types mapped to shapes, and minor
allele frequencies (MAF) mapped to point size.

## Usage

``` r
plot_variant_hotspot(
  var_df,
  chrom_col = "chrom.x",
  pos_col = "pos.x",
  maf_col = "minor_allele_freq",
  start_pos = NULL,
  end_pos = NULL
)
```

## Arguments

- var_df:

  A data frame containing variant data. Required columns: `variant_type`
  (character: "SNP", "INDEL"), `variant_id`, `impact` (character:
  "MODIFIER", "LOW", "MODERATE", "HIGH"), and a numeric MAF column.

- chrom_col:

  A character string specifying the chromosome column.

- pos_col:

  A character string specifying the variant position column.

- maf_col:

  A character string specifying the minor allele frequency column.

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
<<<<<<< HEAD
# Extract genotypes within the candidate gene: Sobic.005G213600
gt_region <- fetch_table_region(table_name = "genotypes",
=======
# Extract variants within the candidate gene: Sobic.005G213600
pg_gt_region <- pg_query_db(table_name = "genotypes",
                            chrom = "Chr05",
                            start = 75104537,
                            end = 75106403
                            )

#' # Extract annotations for the candidate gene
pg_annota_region <- pg_query_db(table_name = "annotations",
>>>>>>> upstream/main
                                chrom = "Chr05",
                                start = 75104537,
                                end = 75106403,
                                connect_db_mode = "online"
                                )

<<<<<<< HEAD
# Extract annotations for the candidate gene
annota_region <- fetch_table_region(table_name = "annotations",
                                    chrom = "Chr05",
                                    start = 75104537,
                                    end = 75106403,
                                    gene_name = "Sobic.005G213600",
                                    connect_db_mode = "online"
                                    )

var_df <- merge(annota_region, gt_region, by = "variant_id")
=======
var_df <- merge(pg_annota_region, pg_gt_region, by = "variant_id")[,1:23]
>>>>>>> upstream/main

plot_variant_hotspot(var_df)

# }
```
