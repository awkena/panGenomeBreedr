# Package index

## Variant Discovery & KASP Marker Design

- [`query_db()`](https://awkena.github.io/panGenomeBreedr/reference/query_db.md)
  : Query any table in your SQLite database using chromosome and a
  genomic position range.
- [`query_by_af()`](https://awkena.github.io/panGenomeBreedr/reference/query_by_af.md)
  : Extract variants based on minimum and maximum allele frequencies
  within a defined region in a SQLite database.
- [`query_ann_summary()`](https://awkena.github.io/panGenomeBreedr/reference/query_ann_summary.md)
  : Query the annotations table within a specified genomic region and
  summarize the distribution of SnpEff annotations and impact categories
  by variant type.
- [`kasp_marker_design()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_marker_design.md)
  : Design KASP markers based on causal variants.
- [`list_sqlite_tables()`](https://awkena.github.io/panGenomeBreedr/reference/list_sqlite_tables.md)
  : List all tables in the SQLite database.
- [`variant_stats()`](https://awkena.github.io/panGenomeBreedr/reference/variant_stats.md)
  : Get variants statistics stored in SQLite database
- [`variant_impact_summary()`](https://awkena.github.io/panGenomeBreedr/reference/variant_impact_summary.md)
  : Get variants statistics stored in SQLite database based on mutation
  impact.
- [`summarize_sqlite_tables()`](https://awkena.github.io/panGenomeBreedr/reference/summarize_sqlite_tables.md)
  : Name and row count for each table in SQLite database.
- [`list_table_columns()`](https://awkena.github.io/panGenomeBreedr/reference/list_table_columns.md)
  : Check the column names and types for any table in a SQLite database.
- [`gene_coord_gff()`](https://awkena.github.io/panGenomeBreedr/reference/gene_coord_gff.md)
  : Get the genomic range of a candidate gene using the Sobic ID from a
  gff file.
- [`query_by_impact()`](https://awkena.github.io/panGenomeBreedr/reference/query_by_impact.md)
  : Extract variants from annotation table based on impact type: LOW,
  MODERATE, HIGH, MODIFIER.
- [`calc_af()`](https://awkena.github.io/panGenomeBreedr/reference/calc_af.md)
  : Compute allele frequencies for a VCF genotype matrix (variant x
  samples). Chromosome and position may be included in the data.
- [`filter_by_af()`](https://awkena.github.io/panGenomeBreedr/reference/filter_by_af.md)
  : Filter extracted variants based on alternate allele frequency.
- [`query_genotypes()`](https://awkena.github.io/panGenomeBreedr/reference/query_genotypes.md)
  : Query genotypes for one or more variant IDs from a wide-format
  genotype table.
- [`count_variant_types()`](https://awkena.github.io/panGenomeBreedr/reference/count_variant_types.md)
  : Count the number of variant types in the SQLite database.
- [`extract_variant()`](https://awkena.github.io/panGenomeBreedr/reference/extract_variant.md)
  : Extract putative causal variants within a candidate gene from a
  tabix-indexed snpEff annotated VCF file.
- [`get_google_id()`](https://awkena.github.io/panGenomeBreedr/reference/get_google_id.md)
  : Get the folder or file ID from a Google Drive shareable link.
- [`folder_download_gd()`](https://awkena.github.io/panGenomeBreedr/reference/folder_download_gd.md)
  : Download files in a shared Google Drive folder without restrictions.

## KASP Marker Validation

- [`read_kasp_csv()`](https://awkena.github.io/panGenomeBreedr/reference/read_kasp_csv.md)
  : Read raw KASP results file (csv format) with one or multiple plates.
- [`get_alleles()`](https://awkena.github.io/panGenomeBreedr/reference/get_alleles.md)
  : Get SNP or InDel alleles and possible genotypes from genotype calls
  in KASP data.
- [`kasp_pch()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_pch.md)
  : Generate pch characters for cluster plots of KASP genotype calls.
- [`kasp_color()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_color.md)
  : Color-code KASP genotype calls based on LGC genomics colors for HEX
  and FAM.
- [`scale_axis()`](https://awkena.github.io/panGenomeBreedr/reference/scale_axis.md)
  : Normalize FAM and HEX fluorescence values between 0 and 1
- [`pred_status()`](https://awkena.github.io/panGenomeBreedr/reference/pred_status.md)
  : Generate the prediction status of positive controls in a KASP assay,
  if present.
- [`pred_summary()`](https://awkena.github.io/panGenomeBreedr/reference/pred_summary.md)
  : Generate summary of prediction for positive controls in KASP
  genotype data, if present
- [`kasp_qc_ggplot()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_qc_ggplot.md)
  : Make KASP marker genotyping QC plot.
- [`kasp_qc_ggplot2()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_qc_ggplot2.md)
  : Make KASP marker genotyping QC plot overlaid with predicitons.
- [`plot_plate()`](https://awkena.github.io/panGenomeBreedr/reference/plot_plate.md)
  : Plot kasp genotyping plate layout.
- [`nsamples_plate()`](https://awkena.github.io/panGenomeBreedr/reference/nsamples_plate.md)
  : Get a summary of the number of samples per 96-well plate in a
  multi-plate KASP assay.
- [`kasp_reshape_wide()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_reshape_wide.md)
  : Reshape KASP data to wide format for same samples genotyped with
  multiple KASP markers.
- [`proc_kasp()`](https://awkena.github.io/panGenomeBreedr/reference/proc_kasp.md)
  : Process reshaped KASP genotype data for heatmap plotting
- [`geno_error()`](https://awkena.github.io/panGenomeBreedr/reference/geno_error.md)
  : Identify SNP loci with potential genotype call errors.
- [`kasp_numeric()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_numeric.md)
  : Convert processed KASP data to numeric genotypes
- [`pred_summary_plot()`](https://awkena.github.io/panGenomeBreedr/reference/pred_summary_plot.md)
  : Create decision support bar plots of match vs. mismatch rates of
  KASP markers that had predictions for positive controls.
- [`gg_dat()`](https://awkena.github.io/panGenomeBreedr/reference/gg_dat.md)
  : Convert data matrix for genotypes to a long format data frame.

## Decision Support for MABC and Trait Introgression

- [`parse_marker_ns()`](https://awkena.github.io/panGenomeBreedr/reference/parse_marker_ns.md)
  : Parse marker names of Hapmap format with a common pattern containing
  chromosome numbers and positions into a map file.
- [`cross_qc_ggplot()`](https://awkena.github.io/panGenomeBreedr/reference/cross_qc_ggplot.md)
  : Create a heatmap to visualize and compare the genetic genetic
  backgrounds of genotypes/lines.
- [`rm_mono()`](https://awkena.github.io/panGenomeBreedr/reference/rm_mono.md)
  : Remove or filter out monomorphic loci from a data matrix or frame.
- [`calc_rpp_bc()`](https://awkena.github.io/panGenomeBreedr/reference/calc_rpp_bc.md)
  : Calculate the proportion of recurrent parent background (RPP) fully
  recovered in backcross progenies.
- [`calc_rpp_exp()`](https://awkena.github.io/panGenomeBreedr/reference/calc_rpp_exp.md)
  : Compute theoretical RPP values for any specified backcross
  generation.
- [`rpp_barplot()`](https://awkena.github.io/panGenomeBreedr/reference/rpp_barplot.md)
  : Visualize computed RPP values for BC progenies as a bar plot.
- [`cross_qc_annotate()`](https://awkena.github.io/panGenomeBreedr/reference/cross_qc_annotate.md)
  : Annotate start and end positions of loci on a heatmap.
- [`cross_qc_heatmap()`](https://awkena.github.io/panGenomeBreedr/reference/cross_qc_heatmap.md)
  : Create a heatmap to visualize and compare the genetic genetic
  backgrounds of genotypes/lines with or without annotation for
  introgressed loci.
- [`cross_qc_heatmap2()`](https://awkena.github.io/panGenomeBreedr/reference/cross_qc_heatmap2.md)
  : Visualize genotype backgrounds with optional QTL annotations.
- [`sim_snp_dat()`](https://awkena.github.io/panGenomeBreedr/reference/sim_snp_dat.md)
  : Simulate raw SNP loci for any chromosome with or without LD.
- [`order_markers()`](https://awkena.github.io/panGenomeBreedr/reference/order_markers.md)
  : Order marker IDs based on their chromosome numbers and positions in
  ascending order.
- [`hapmap_ns_fmt()`](https://awkena.github.io/panGenomeBreedr/reference/hapmap_ns_fmt.md)
  : Format marker names to comply with the Hapmap convention.
- [`find_unexp_homs()`](https://awkena.github.io/panGenomeBreedr/reference/find_unexp_homs.md)
  : Find loci with unexpected homozygous genotype calls for artificial
  heterozygotes.
- [`find_indels()`](https://awkena.github.io/panGenomeBreedr/reference/find_indels.md)
  : Identify and subset InDel markers from a marker panel.
- [`parent_missing()`](https://awkena.github.io/panGenomeBreedr/reference/parent_missing.md)
  : Identify and subset loci with any parent missing genotype.
- [`parent_het()`](https://awkena.github.io/panGenomeBreedr/reference/parent_het.md)
  : Identify and subset loci with any heterozygous parent genotype.
- [`parent_poly()`](https://awkena.github.io/panGenomeBreedr/reference/parent_poly.md)
  : Select polymorphic loci between two parents in a marker panel.
- [`foreground_select()`](https://awkena.github.io/panGenomeBreedr/reference/foreground_select.md)
  : Identify lines that possess favorable alleles for target loci using
  trait predictive markers.
- [`find_lines()`](https://awkena.github.io/panGenomeBreedr/reference/find_lines.md)
  : Extracts lines that have a combination of favorable alleles across
  target loci.

## Example Datasets

- [`beta_carotene`](https://awkena.github.io/panGenomeBreedr/reference/beta_carotene.md)
  : beta_carotene
- [`kasp_dat`](https://awkena.github.io/panGenomeBreedr/reference/kasp_dat.md)
  : kasp_dat

## Applications

- [`run_app()`](https://awkena.github.io/panGenomeBreedr/reference/run_app.md)
  : Run the Shiny Application
