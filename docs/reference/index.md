# Package index

## Variant Discovery & KASP Marker Design

### Connection & Configuration

- [`connect_local_db()`](https://awkena.github.io/panGenomeBreedr/reference/connect_local_db.md)
  : Connect to the local offline database
- [`disconnect_local_db()`](https://awkena.github.io/panGenomeBreedr/reference/disconnect_local_db.md)
  : Disconnect the local database session
- [`set_api_url()`](https://awkena.github.io/panGenomeBreedr/reference/set_api_url.md)
  : Set the API endpoint URL for your private database
- [`get_api_url()`](https://awkena.github.io/panGenomeBreedr/reference/get_api_url.md)
  : Get the current database API endpoint URL

### Database Schema & Metadata

- [`list_tables()`](https://awkena.github.io/panGenomeBreedr/reference/list_tables.md)
  : List all tables in the database
- [`summarize_database()`](https://awkena.github.io/panGenomeBreedr/reference/summarize_database.md)
  : Get table names and row counts from the database
- [`list_columns()`](https://awkena.github.io/panGenomeBreedr/reference/list_columns.md)
  : Get table schema details

### Core Data Queries

- [`fetch_table_region()`](https://awkena.github.io/panGenomeBreedr/reference/fetch_table_region.md)
  : Fetch data from database tables by genomic coordinates
- [`fetch_genotypes_by_id()`](https://awkena.github.io/panGenomeBreedr/reference/fetch_genotypes_by_id.md)
  : Fetch genotypes for specific variant IDs
- [`fetch_accession_metadata()`](https://awkena.github.io/panGenomeBreedr/reference/fetch_accession_metadata.md)
  : Get sample metadata from the database
- [`filter_genotypes_by_metadata()`](https://awkena.github.io/panGenomeBreedr/reference/filter_genotypes_by_metadata.md)
  : Filter a genotype matrix by sample metadata

### Targeted Variant Filters

- [`fetch_variants_by_allele_frequency()`](https://awkena.github.io/panGenomeBreedr/reference/fetch_variants_by_allele_frequency.md)
  : Fetch variants by allele frequency
- [`fetch_variants_by_impact()`](https://awkena.github.io/panGenomeBreedr/reference/fetch_variants_by_impact.md)
  : Extract variants based on functional mutation impact

### Summaries & Statistics

- [`summarize_variants()`](https://awkena.github.io/panGenomeBreedr/reference/summarize_variants.md)
  : Get variant statistics from the database
- [`summarize_variant_impacts()`](https://awkena.github.io/panGenomeBreedr/reference/summarize_variant_impacts.md)
  : Get summary of database variant impacts by chromosome
- [`count_variant_types()`](https://awkena.github.io/panGenomeBreedr/reference/count_variant_types.md)
  : Count the distribution of variant types in the database
- [`summarize_annotations()`](https://awkena.github.io/panGenomeBreedr/reference/summarize_annotations.md)
  : Summarize functional annotations and impacts for a genomic region

### Allele Frequency Handlers

- [`calculate_allele_frequencies()`](https://awkena.github.io/panGenomeBreedr/reference/calculate_allele_frequencies.md)
  : Compute allele frequencies for a genotype matrix
- [`filter_by_allele_frequency()`](https://awkena.github.io/panGenomeBreedr/reference/filter_by_allele_frequency.md)
  : Filter genotypes by allele frequency

### Linkage Disequilibrium Analysis

- [`calculate_LD()`](https://awkena.github.io/panGenomeBreedr/reference/calculate_LD.md)
  : Compute Linkage Disequilibrium (LD) Metrics (R2 and D')
- [`plot_ld_geodesic()`](https://awkena.github.io/panGenomeBreedr/reference/plot_ld_geodesic.md)
  : Generate Geodesic Landscape Plot and Extract Haploblock Data

### Visualization

- [`plot_accession_map()`](https://awkena.github.io/panGenomeBreedr/reference/plot_accession_map.md)
  : Interactive geographic exploration of the samples
- [`plot_variant_hotspot()`](https://awkena.github.io/panGenomeBreedr/reference/plot_variant_hotspot.md)
  : Plot genomic variant hotspots (SNPs and INDELs) extracted from the
  database (local)
- [`plot_gene_model()`](https://awkena.github.io/panGenomeBreedr/reference/plot_gene_model.md)
  : Plot the gene model using genomic coordinates (local)
- [`hotspot_overlay_plot()`](https://awkena.github.io/panGenomeBreedr/reference/hotspot_overlay_plot.md)
  : Generate a combined gene model and variant hotspot overlay plot
  (local)

### Marker Design & Utilities

- [`kasp_marker_design()`](https://awkena.github.io/panGenomeBreedr/reference/kasp_marker_design.md)
  : Design KASP markers based on causal variants.
- [`gene_coord_gff()`](https://awkena.github.io/panGenomeBreedr/reference/gene_coord_gff.md)
  : Get the genomic range of a candidate gene and its features using the
  Sobic ID from a GFF file (local)
- [`extract_variant()`](https://awkena.github.io/panGenomeBreedr/reference/extract_variant.md)
  : Extract putative causal variants within a candidate gene from a
  tabix-indexed snpEff annotated VCF file.
- [`folder_download_gd()`](https://awkena.github.io/panGenomeBreedr/reference/folder_download_gd.md)
  : Download files in a shared Google Drive folder without restrictions.
- [`get_google_id()`](https://awkena.github.io/panGenomeBreedr/reference/get_google_id.md)
  : Get the folder or file ID from a Google Drive shareable link.

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
