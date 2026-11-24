#' Reproducible-script steps for the Variant Discovery tab
#'
#' Each function takes the settings the tab used, as a plain list, and
#' returns a `repro_step()`. Unlike the other tabs, the script does not come
#' with a data bundle: it queries the same database the app was connected to
#' (a local Parquet folder or the online API), so it needs access to it.
#'
#' Objects created by the script, for the steps that build on each other:
#' `con` (local database connection), `region` (chrom, start, end),
#' `variants`, `annotations`, `genotypes` (the region's tables),
#' `genotypes_shown`, `variants_shown`, `annotations_shown` (after the
#' table filters), `selected_variants` (variants kept by "Refine plot", or
#' NULL), `pcil_data`, `pcil_families`, `pcil_positive`.
#'
#' @name repro_vd
#' @keywords internal
#' @noRd
NULL


#' Header title of the Variant Discovery script
#' @noRd
repro_vd_title <- "panGenomeBreedr variant discovery"


#' Script file name for the Variant Discovery tab
#' @noRd
repro_vd_script <- "variant_discovery.R"


#' Header lines on where the script gets its data
#' @noRd
repro_vd_notes <- c(
  "This script queries the same database as the app (the local folder or",
  "online API set in step 1), so it needs access to it. Results are",
  "written to results/."
)


#' Plain names of the steps, shown when a step still has to run
#' @noRd
repro_vd_labels <- c(
  connect = "Connect to a database",
  region = "Identify Variants: define a genomic region",
  filters = "Identify Variants: load a region",
  pcil_data = "PCIL Selection: load the PCIL data",
  pcil_families = "PCIL Selection: find donor families again",
  pcil_positive = "PCIL Selection: identify the positive lines again",
  pcil_negative = "PCIL Selection: identify the negative lines"
)


#' Columns of the genotype table that are not samples
#' @noRd
repro_vd_info_cols <- c(
  "variant_id", "chrom", "pos", "ref", "alt", "variant_type",
  "major_allele", "minor_allele", "major_allele_freq", "minor_allele_freq"
)


#' Database argument of the panGB query functions
#'
#' @param mode `"local"` or `"online"`.
#' @return A list to put first in a `repro_call()` argument list.
#' @noRd
repro_vd_db <- function(mode) {
  if (identical(mode, "online")) {
    list(connect_db_mode = "online")
  } else {
    list(con = repro_sym("con"))
  }
}


#' Step: connect to the database
#'
#' @param s List with `mode` (`"local"` or `"online"`), and `db_path` and
#'   `example` (TRUE for the bundled example database) for a local
#'   database, or `api_url` for the online one.
#' @noRd
repro_vd_connect_step <- function(s) {
  if (identical(s$mode, "online")) {
    code <- c(
      "# Online database (API) used in the app",
      repro_call("set_api_url", list(url = s$api_url))
    )
    settings <- list(mode = "online", api_url = s$api_url)
  } else if (isTRUE(s$example)) {
    code <- c(
      "# The example database that comes with panGenomeBreedr",
      "db_path <- system.file(\"extdata\", \"pangenome_scale_db\", package = \"panGenomeBreedr\")",
      "con <- connect_local_db(folder_path = db_path)"
    )
    settings <- list(mode = "local", database = "example (pangenome_scale_db)")
  } else {
    code <- c(
      "# Local database folder used in the app. Change the path if you run",
      "# this on another computer.",
      paste("db_path <-", repro_value(s$db_path)),
      "con <- connect_local_db(folder_path = db_path)"
    )
    settings <- list(mode = "local", db_path = s$db_path)
  }

  repro_step(
    title = "Connect to the database",
    code = c(
      code,
      "",
      "# Save a plot, or a list of plots, to one PDF in results/",
      "save_plots <- function(plots, file, width = 12, height = 8) {",
      "  if (inherits(plots, \"ggplot\")) plots <- list(plots)",
      "  pdf(file.path(\"results\", file), width = width, height = height, onefile = TRUE)",
      "  for (p in plots) print(p)",
      "  invisible(dev.off())",
      "}"
    ),
    settings = settings
  )
}


#' Step: close the local database connection
#' @noRd
repro_vd_disconnect_step <- function() {
  repro_step(
    title = "Close the database connection",
    code = "disconnect_local_db(con)",
    needs = "connect"
  )
}


#' Code that sets `gff_file`, the GFF3 source of a gene search
#' @noRd
repro_vd_gff_code <- function(s) {
  if (identical(s$gff_source, "file")) {
    name <- repro_file_name(if (is.null(s$gff_name)) "genes.gff3" else s$gff_name)
    c(
      "# The GFF3 file you uploaded in the app. Put it next to this script,",
      "# or change the path.",
      paste("gff_file <-", repro_value(name)),
      "if (!file.exists(gff_file)) stop(\"Can't find the GFF3 file \", gff_file)"
    )
  } else {
    paste("gff_file <-", repro_value(s$gff_url))
  }
}


#' Step: define the region and fetch its tables
#'
#' @param s List with `method` (`"gene"` or `"coords"`), `mode`, and either
#'   `gene_name`, `gff_source` (`"url"` or `"file"`), `gff_url`,
#'   `gff_name`, or `chrom`, `start`, `end`.
#' @noRd
repro_vd_region_step <- function(s) {
  db <- repro_vd_db(s$mode)
  reg <- list(
    chrom = repro_sym("region$chrom"),
    start = repro_sym("region$start"),
    end = repro_sym("region$end")
  )

  if (identical(s$method, "gene")) {
    define <- c(
      repro_vd_gff_code(s),
      repro_call("gene_coord_gff",
                 list(gene_name = s$gene_name, gff_path = repro_sym("gff_file")),
                 assign = "gene_models"),
      "gene_rows <- gene_models[gene_models$Feature == \"gene\", ]",
      "region <- list(",
      "  chrom = gene_rows$Chromosome[1],",
      "  start = min(gene_rows$Start, na.rm = TRUE),",
      "  end = max(gene_rows$End, na.rm = TRUE)",
      ")"
    )
    settings <- list(method = "gene", gene_name = s$gene_name,
                     gff = if (identical(s$gff_source, "file")) s$gff_name else s$gff_url)
  } else {
    define <- paste("region <-", repro_value(list(chrom = s$chrom, start = s$start, end = s$end)))
    settings <- list(method = "coordinates", chrom = s$chrom, start = s$start, end = s$end)
  }

  fetch <- unlist(lapply(c("variants", "annotations", "genotypes"), function(tb) {
    repro_call("fetch_table_region", c(db, list(table_name = tb), reg), assign = tb)
  }))

  repro_step(
    title = "Define the genomic region and fetch its variants",
    code = c(
      define,
      "",
      fetch,
      repro_call("summarize_annotations", c(db, reg), assign = "annotation_summary"),
      "",
      repro_save_table("variants", "region_variants.csv"),
      repro_save_table("annotations", "region_annotations.csv"),
      repro_save_table("genotypes", "region_genotypes.csv"),
      repro_save_table("annotation_summary$annotation_summary", "annotation_summary.csv"),
      repro_save_table("annotation_summary$impact_summary", "impact_summary.csv"),
      repro_save_table("annotation_summary$variant_type_totals", "variant_type_totals.csv")
    ),
    settings = settings,
    needs = "connect"
  )
}


#' A filter list, written as R code with one filter per line
#' @noRd
repro_vd_filter_code <- function(filters, name) {
  tables <- c("genotypes", "variants", "annotations")
  lines <- vapply(tables, function(tb) {
    f <- filters[[tb]]
    if (length(f) == 0) f <- list()
    f <- lapply(f, function(x) list(column = x$column, value = x$value))
    paste0("  ", tb, " = ", repro_value(f))
  }, character(1))
  c(paste(name, "<- list("),
    paste0(lines, c(",", ",", "")),
    ")")
}


#' TRUE if a list of table filters holds at least one filter
#' @noRd
repro_vd_has_filters <- function(filters) {
  any(lengths(filters[c("genotypes", "variants", "annotations")]) > 0)
}


#' Step: table filters and the variants kept by "Refine plot"
#'
#' @param s List with `mode`, `filters` (list of `genotypes`, `variants`,
#'   `annotations`, each a list of `list(column, value)`), `sample_filters`
#'   (named list: metadata column = values to keep) and `plot_filters`
#'   (the filters when "Refine plot" was last clicked, or NULL).
#' @noRd
repro_vd_filters_step <- function(s) {
  filters <- s$filters
  plot_filters <- s$plot_filters
  has_filters <- repro_vd_has_filters(filters)
  has_samples <- length(s$sample_filters) > 0
  has_plot <- !is.null(plot_filters) && repro_vd_has_filters(plot_filters)

  if (!has_filters && !has_samples && !has_plot) {
    return(repro_step(
      title = "Table filters",
      code = c(
        "# No filters were applied in the app: keep every variant and sample",
        "genotypes_shown <- genotypes",
        "variants_shown <- variants",
        "annotations_shown <- annotations",
        "selected_variants <- NULL"
      ),
      settings = list(filters = "none"),
      needs = "region"
    ))
  }

  code <- c(
    "# Keep the rows that pass each filter: a (min, max) range for numeric",
    "# columns, or the values to keep for the others",
    "apply_filters <- function(df, filters) {",
    "  for (f in filters) {",
    "    if (!f$column %in% names(df)) next",
    "    df <- df[!is.na(df[[f$column]]), ]",
    "    if (is.numeric(f$value) && is.numeric(df[[f$column]])) {",
    "      df <- df[df[[f$column]] >= f$value[1] & df[[f$column]] <= f$value[2], ]",
    "    } else {",
    "      df <- df[as.character(df[[f$column]]) %in% as.character(f$value), ]",
    "    }",
    "  }",
    "  df",
    "}",
    "",
    repro_vd_filter_code(filters, "filters"),
    "",
    "genotypes_shown <- apply_filters(genotypes, filters$genotypes)",
    "variants_shown <- apply_filters(variants, filters$variants)",
    "annotations_shown <- apply_filters(annotations, filters$annotations)"
  )

  if (has_samples) {
    code <- c(
      code,
      "",
      "# Samples kept by the sample (column) filters",
      paste("sample_filters <-", repro_value(s$sample_filters)),
      repro_call("fetch_accession_metadata", repro_vd_db(s$mode), assign = "samples"),
      "for (col in names(sample_filters)) {",
      "  samples <- samples[samples[[col]] %in% sample_filters[[col]], , drop = FALSE]",
      "}",
      paste("info_cols <-", repro_value(repro_vd_info_cols)),
      "genotypes_shown <- genotypes_shown[, c(",
      "  intersect(info_cols, names(genotypes_shown)),",
      "  intersect(samples$lib, names(genotypes_shown))",
      "), drop = FALSE]"
    )
  }

  if (has_plot) {
    same <- isTRUE(all.equal(plot_filters[c("genotypes", "variants", "annotations")],
                             filters[c("genotypes", "variants", "annotations")]))
    plot_code <- if (same) {
      "plot_filters <- filters"
    } else {
      c("# Filters in place when \"Refine plot\" was last clicked",
        repro_vd_filter_code(plot_filters, "plot_filters"))
    }
    code <- c(
      code,
      "",
      "# Variants left in every filtered table (\"Refine plot\")",
      plot_code,
      "tables <- list(genotypes = genotypes, variants = variants, annotations = annotations)",
      "id_sets <- list()",
      "for (tb in names(tables)) {",
      "  if (length(plot_filters[[tb]]) > 0) {",
      "    id_sets[[tb]] <- apply_filters(tables[[tb]], plot_filters[[tb]])$variant_id",
      "  }",
      "}",
      "selected_variants <- Reduce(intersect, id_sets)"
    )
  } else {
    code <- c(code, "", "# \"Refine plot\" was not used: the plot shows every variant",
              "selected_variants <- NULL")
  }

  code <- c(
    code,
    "",
    repro_save_table("genotypes_shown", "filtered_genotypes.csv"),
    repro_save_table("variants_shown", "filtered_variants.csv"),
    repro_save_table("annotations_shown", "filtered_annotations.csv")
  )

  repro_step(
    title = "Table filters",
    code = code,
    settings = list(
      filters = filters,
      sample_filters = if (has_samples) s$sample_filters,
      plot_filters = if (has_plot) plot_filters
    ),
    needs = "region"
  )
}


#' Step: variant hotspot plot of the region
#'
#' @param s List with `method` (`"gene"` or `"coords"`), `gene_name`,
#'   `text_sz` and `plot_annotations` (TRUE if "Refine plot" filtered the
#'   annotations).
#' @noRd
repro_vd_hotspot_step <- function(s) {
  code <- c(
    "plot_ann <- annotations",
    "plot_geno <- genotypes",
    "if (!is.null(selected_variants)) {",
    "  plot_ann <- plot_ann[plot_ann$variant_id %in% selected_variants, ]",
    "  plot_geno <- plot_geno[plot_geno$variant_id %in% selected_variants, ]",
    "}"
  )
  if (isTRUE(s$plot_annotations)) {
    code <- c(code, "plot_ann <- apply_filters(plot_ann, plot_filters$annotations)")
  }
  code <- c(
    code,
    "",
    "# Keep the most severe impact of each variant",
    "if (nrow(plot_ann) > 0 && \"impact\" %in% names(plot_ann)) {",
    "  impact_rank <- factor(plot_ann$impact,",
    "                        levels = c(\"HIGH\", \"MODERATE\", \"LOW\", \"MODIFIER\"),",
    "                        ordered = TRUE)",
    "  plot_ann <- plot_ann[order(impact_rank, na.last = TRUE), ]",
    "  plot_ann <- plot_ann[!duplicated(plot_ann$variant_id), ]",
    "}",
    ""
  )

  plot_call <- if (identical(s$method, "gene")) {
    c(
      "# hotspot_overlay_plot() is the app's gene-model overlay (internal to panGB)",
      repro_call("panGenomeBreedr:::hotspot_overlay_plot", list(
        gene_name = s$gene_name,
        gff_path = repro_sym("gff_file"),
        annotations_df = repro_sym("plot_ann"),
        genotypes_df = repro_sym("plot_geno"),
        selected_variants = repro_sym("unique(selected_variants)"),
        text_sz = s$text_sz
      ), assign = "hotspot")
    )
  } else {
    repro_call("plot_variant_hotspot", list(
      var_df = repro_sym("merge(plot_ann, plot_geno, by = \"variant_id\")"),
      start_pos = repro_sym("region$start"),
      end_pos = repro_sym("region$end")
    ), assign = "hotspot")
  }

  repro_step(
    title = "Variant hotspot plot",
    code = c(code, plot_call, "save_plots(hotspot, \"hotspot_plot.pdf\", width = 12, height = 8)"),
    settings = list(method = s$method, gene_label_size = s$text_sz),
    needs = "filters"
  )
}


#' Step: linkage disequilibrium analysis
#' @param s List with `target_ids`, `metric`, `threshold`,
#'   `block_threshold`, `filled` and `palette_option`.
#' @noRd
repro_vd_ld_step <- function(s) {
  targets <- if (length(s$target_ids)) s$target_ids else NULL
  repro_step(
    title = "Linkage disequilibrium (LD) analysis",
    code = c(
      repro_call("calculate_LD", list(
        df = repro_sym("genotypes_shown"),
        target_variant_ids = targets
      ), assign = "ld_res"),
      repro_call("plot_ld_geodesic", list(
        ld_df = repro_sym("ld_res"),
        query_db_geno = repro_sym("genotypes_shown"),
        query_db_annot = repro_sym("annotations"),
        metric = s$metric,
        target_variant_ids = targets,
        threshold = s$threshold,
        block_threshold = s$block_threshold,
        show_variant_labels = TRUE,
        filled = s$filled,
        palette_option = s$palette_option
      ), assign = "ld"),
      "save_plots(ld$plot, \"ld_geodesic_plot.pdf\", width = 14, height = 9)",
      repro_save_table("ld$table", "haploblocks.csv"),
      repro_save_table("ld_res", "ld_pairs.csv")
    ),
    settings = list(target_variants = targets, metric = s$metric,
                    threshold = s$threshold, block_threshold = s$block_threshold,
                    filled = s$filled, palette = s$palette_option),
    needs = "filters"
  )
}


#' Step: load the PCIL data
#' @param s List with `mode`.
#' @noRd
repro_vd_pcil_data_step <- function(s) {
  repro_step(
    title = "Load the PCIL data",
    code = repro_call("fetch_pcil_data", repro_vd_db(s$mode), assign = "pcil_data"),
    needs = "connect"
  )
}


#' Step: PCIL donor families carrying the selected variants
#' @param s List with `mode` and `selection` (variant IDs).
#' @noRd
repro_vd_pcil_families_step <- function(s) {
  db <- if (identical(s$mode, "online")) list() else list(con = repro_sym("con"))
  repro_step(
    title = "PCIL donor families",
    code = c(
      repro_call("fetch_pcil_families_by_variant", c(db, list(
        selection = s$selection,
        pcil_data = repro_sym("pcil_data"),
        connect_db_mode = if (identical(s$mode, "online")) "online" else "local"
      )), assign = "pcil_families"),
      repro_save_table("pcil_families$pcil_family_summary", "pcil_family_summary.csv"),
      repro_save_table("pcil_families$pcil_summary", "pcil_families.csv")
    ),
    settings = list(variants = s$selection),
    needs = "pcil_data"
  )
}


#' Code for the `available_ids` and `global_available_ids` of the PCIL calls
#' @noRd
repro_vd_pcil_ids <- function(s) {
  list(
    available_ids = if (isTRUE(s$use_family_ids)) {
      repro_sym("pcil_families$pcil_summary[, c(\"sample_id\", \"selection\")]")
    },
    result_pcil_families = if (isTRUE(s$use_families)) repro_sym("pcil_families"),
    global_available_ids = if (length(s$global_ids)) s$global_ids
  )
}


#' Step: positive PCIL lines
#'
#' @param s List with `type` (`"gene"`, `"region"` or `"position"`),
#'   `gene_ids`, `region` (list of chr, start, end), `variant_ids`, `sel`,
#'   `donor_thresh`, `block_quantile`, `F_quantile`, `window`,
#'   `use_families`, `use_family_ids` and `global_ids`.
#' @noRd
repro_vd_pcil_positive_step <- function(s) {
  needs <- "pcil_data"
  if (identical(s$type, "gene")) {
    targets <- paste("pcil_targets <-", repro_value(s$gene_ids))
    target_settings <- list(genes = s$gene_ids)
  } else if (identical(s$type, "region")) {
    r <- s$region
    region_id <- sprintf("%s:%s-%s", r$chr,
                         format(r$start, scientific = FALSE, trim = TRUE),
                         format(r$end, scientific = FALSE, trim = TRUE))
    targets <- c(
      "pcil_targets <- data.frame(",
      paste0("  Region = ", repro_value(region_id), ","),
      paste0("  Chr = ", repro_value(r$chr), ","),
      paste0("  Start = ", repro_value(r$start), ","),
      paste0("  End = ", repro_value(r$end), ","),
      "  stringsAsFactors = FALSE",
      ")"
    )
    target_settings <- list(region = r)
  } else {
    targets <- c(
      "# The region's variants chosen in the app",
      paste0("pcil_targets <- variants[variants$variant_id %in% ",
             repro_value(s$variant_ids), ", ]")
    )
    target_settings <- list(variants = s$variant_ids)
    needs <- c(needs, "region")
  }
  if (isTRUE(s$use_families)) needs <- c(needs, "pcil_families")

  args <- c(
    list(
      pcil_data = repro_sym("pcil_data"),
      variants_select_geno = repro_sym("pcil_targets"),
      type = s$type,
      sel = s$sel,
      donor_thresh = s$donor_thresh,
      block_quantile = s$block_quantile,
      F_quantile = s$F_quantile,
      window = s$window
    ),
    repro_vd_pcil_ids(s)
  )

  repro_step(
    title = "Positive PCIL lines",
    code = c(
      targets,
      repro_call("fetch_pcil_positive", args, assign = "pcil_positive"),
      repro_save_table("pcil_positive$best_lines", "pcil_positive_best_lines.csv"),
      repro_save_table("pcil_positive$pcil_positive", "pcil_positive_lines.csv"),
      "",
      "save_plots(plot_pcil_positive(pcil_positive_result = pcil_positive),",
      "           \"pcil_positive_per_region.pdf\")",
      "save_plots(plot_pcil_best_lines(pcil_positive_result = pcil_positive, pcil_data = pcil_data),",
      "           \"pcil_positive_best_lines.pdf\")"
    ),
    settings = c(list(type = s$type), target_settings, list(
      sel = s$sel, donor_thresh = s$donor_thresh,
      block_quantile = s$block_quantile, F_quantile = s$F_quantile,
      window = s$window, available_ids = if (length(s$global_ids)) s$global_ids
    )),
    needs = needs
  )
}


#' Step: negative control PCIL lines
#' @param s List with `n_neg`, `use_families`, `use_family_ids` and
#'   `global_ids`.
#' @noRd
repro_vd_pcil_negative_step <- function(s) {
  args <- c(
    list(
      pcil_data = repro_sym("pcil_data"),
      pcil_positive_result = repro_sym("pcil_positive"),
      n_neg = s$n_neg
    ),
    repro_vd_pcil_ids(s)
  )
  repro_step(
    title = "Negative control PCIL lines",
    code = c(
      repro_call("fetch_pcil_negative", args, assign = "pcil_negative"),
      repro_save_table("pcil_negative$pairs_best", "pcil_negative_best.csv"),
      repro_save_table("pcil_negative$pairs_extended", "pcil_negative_ranked.csv"),
      "save_plots(plot_pcil_pairs(pcil_neg_results = pcil_negative, pcil_data = pcil_data),",
      "           \"pcil_negative_pairs.pdf\")"
    ),
    settings = list(n_neg = s$n_neg),
    needs = c("pcil_positive", if (isTRUE(s$use_families)) "pcil_families")
  )
}


#' Step: KASP markers for the region's variants
#'
#' @param s List with `mode`, `candidates` (variant IDs sent to marker
#'   design, or NULL for every variant in the region), `meta_cols`,
#'   `genome_name`, `design` (arguments of kasp_marker_design() other than
#'   `gt_df`, `chr`, `genome_file`, `plot_file` and `save_alignment`) and
#'   `intertek` (arguments of make_intertek_table() other than
#'   `marker_data`).
#' @noRd
repro_vd_kasp_step <- function(s) {
  genome <- repro_file_name(if (is.null(s$genome_name)) "genome.fa.gz" else s$genome_name)
  candidates <- if (is.null(s$candidates)) {
    "marker_candidates <- unique(variants$variant_id)"
  } else {
    paste("marker_candidates <-", repro_value(s$candidates))
  }
  d <- s$design

  repro_step(
    title = "KASP marker design",
    code = c(
      "# Variants sent to marker design",
      candidates,
      repro_call("fetch_genotypes_by_id", c(repro_vd_db(s$mode), list(
        variant_ids = repro_sym("marker_candidates"),
        meta_data = s$meta_cols
      )), assign = "marker_genotypes"),
      "",
      "# The genome FASTA you uploaded in the app. Put it next to this script,",
      "# or change the path.",
      paste("genome_file <-", repro_value(genome)),
      "if (!file.exists(genome_file)) stop(\"Can't find the genome file \", genome_file)",
      "",
      repro_call("kasp_marker_design", c(
        list(gt_df = repro_sym("marker_genotypes")),
        d[c("variant_id_col", "chrom_col", "pos_col", "ref_al_col", "alt_al_col",
            "geno_start", "marker_IDs")],
        list(
          chr = repro_sym(sprintf("unique(marker_genotypes[[%s]])", repro_value(d$chrom_col))),
          genome_file = repro_sym("genome_file"),
          plot_file = repro_sym("tempdir()")
        ),
        d[c("region_name", "maf")],
        list(save_alignment = FALSE)
      ), assign = "markers"),
      repro_call("make_intertek_table",
                 c(list(marker_data = repro_sym("markers")), s$intertek),
                 assign = "intertek"),
      repro_save_table("intertek", "KASP_markers_for_Intertek.csv"),
      "save_plots(markers$plot, \"alignment_plots.pdf\", width = 24, height = 9)"
    ),
    settings = c(d, s$intertek, list(genome_file = genome)),
    needs = if (is.null(s$candidates)) c("connect", "region") else "connect"
  )
}
