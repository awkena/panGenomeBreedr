#' Reproducible-script steps for the Breeding Decision Support tabs
#'
#' Each function takes the settings a module used, as a plain list, and
#' returns a `repro_step()` (see `R/utils_repro.R`). The scripts spell out
#' the exported panGB calls behind the app's internal helpers, so they can
#' be read and changed without the app.
#'
#' Objects created by the processing step, used by the later steps:
#' `geno` (genotype file), `geno_qc` (after the marker QC filters),
#' `geno_num` (numeric-coded genotypes) and `map` (map file).
#'
#' @name repro_ds
#' @keywords internal
#' @noRd
NULL


#' Script titles and file names, one per tab
#' @noRd
repro_ds_titles <- c(
  introgression = "panGenomeBreedr trait introgression check",
  mabc = "panGenomeBreedr marker-assisted backcrossing (RPP)",
  foreground = "panGenomeBreedr foreground selection"
)

#' @noRd
repro_ds_scripts <- c(
  introgression = "trait_introgression.R",
  mabc = "marker_assisted_backcrossing.R",
  foreground = "foreground_selection.R"
)


#' Code that reads a csv or Excel file the way read_mapfile() does
#'
#' @param obj Object name to assign to.
#' @param name File name in `data/`.
#' @noRd
repro_read_table_code <- function(obj, name) {
  path <- repro_value(repro_data_path(name))
  if (grepl("\\.xlsx?$", name, ignore.case = TRUE)) {
    sprintf("%s <- as.data.frame(readxl::read_excel(%s, .name_repair = \"minimal\"))", obj, path)
  } else {
    sprintf("%s <- read.csv(%s, stringsAsFactors = FALSE)", obj, path)
  }
}


#' Add text to the last line of some code, e.g. `$par_present` after a call
#' @noRd
repro_suffix <- function(code, suffix) {
  code[length(code)] <- paste0(code[length(code)], suffix)
  code
}


#' Write a data frame to a temporary csv, for inputs that never were a file
#' (such as the built-in example of the foreground tab)
#' @noRd
repro_temp_csv <- function(df, name) {
  path <- file.path(tempdir(), repro_file_name(name))
  utils::write.csv(df, path, row.names = FALSE)
  path
}


#' Settings of the Process Data step, read from the module inputs
#'
#' Shared by the trait introgression and MABC modules, which have the same
#' Process Data sidebar. Call it inside a reactive context, right after
#' proc_nd_map_func() succeeds.
#'
#' @param input The module's `input`.
#' @param data_path,data_name Local path and original name of the genotype
#'   file.
#' @param result Value of proc_nd_map_func().
#' @param map_path,map_name Local path and original name of the map file,
#'   if one was uploaded.
#' @param map_data The uploaded map file, read in, if any.
#' @param data_cols Column names of the genotype file.
#' @noRd
repro_ds_process_settings <- function(input, data_path, data_name, result,
                                      map_path = NULL, map_name = NULL,
                                      map_data = NULL, data_cols = NULL) {
  map_cols <- colnames(result$mapfile)
  pick <- function(pattern, cols) safe_grep_match(pattern, choices = cols)[1]
  batch <- !is.null(input$cluster_by) && input$cluster_by != "None" &&
    input$cluster_by %in% data_cols
  list(
    data_name = if (is.null(data_name)) "genotypes.csv" else data_name,
    data_path = data_path,
    cluster_by = if (batch) input$cluster_by,
    cluster_focus = if (batch) input$cluster_focus,
    marker_start = input$marker_start,
    sample_id = input$sample_id,
    calls_sep = check_sep(input$allele_sep),
    data_type = if (input$data_type == "Agriplex") "agriplex" else if (input$data_type == "Kasp / DArTag") "kasp",
    rp = as.numeric(input$rp),
    dp = as.numeric(input$dp),
    rp_name = result$rp_name,
    dp_name = result$dp_name,
    apply_par_poly = isTRUE(input$apply_par_poly),
    apply_par_miss = isTRUE(input$apply_par_miss),
    apply_geno_good = isTRUE(input$apply_geno_good),
    apply_par_homo = isTRUE(input$apply_par_homo),
    feedback = input$choice,
    marker_sep = input$sep_marker,
    prefix = input$prefix_marker,
    map_name = map_name,
    map_path = map_path,
    map_snp_col = if (!is.null(map_data)) pick("snp", colnames(map_data)),
    map_cols = c(snpid = pick("snpid", map_cols), chr = pick("chr", map_cols),
                 pos = pick("pos", map_cols))
  )
}


#' Step: read the genotype file, run the marker QC filters, build the map
#' and code the genotypes as numbers
#'
#' Spells out what the app's Process Data button runs.
#'
#' @param s List with `data_name`, `data_path`, `cluster_by`,
#'   `cluster_focus` (NULL for no batch filter), `marker_start`,
#'   `sample_id`, `calls_sep`, `data_type`, `rp`, `dp` (row numbers),
#'   `rp_name`, `dp_name`, the four `apply_*` filter switches, `feedback`
#'   (`"no"`: map from marker names, `"yes"`: map file), `marker_sep`,
#'   `prefix`, `map_name`, `map_path`, `map_snp_col` and `map_cols`
#'   (named `snpid`, `chr`, `pos`: map columns used by proc_kasp()).
#' @noRd
repro_ds_process_step <- function(s) {
  data_name <- repro_file_name(s$data_name)
  ms <- s$marker_start
  qc_args <- list(x = repro_sym("calls"), rp_row = repro_sym("rp_row"),
                  dp_row = repro_sym("dp_row"))
  sep_arg <- list(sep = repro_sym("calls_sep"))
  any_filter <- any(s$apply_par_poly, s$apply_par_miss, s$apply_geno_good,
                    s$apply_par_homo)

  code <- c(
    "# Genotype calls",
    repro_read_table_code("geno", data_name),
    if (!is.null(s$cluster_focus)) {
      sprintf("geno <- geno[geno[[%s]] == %s, ]  # keep only %s %s",
              repro_value(s$cluster_by), repro_value(as.character(s$cluster_focus)),
              s$cluster_by, s$cluster_focus)
    },
    "",
    sprintf("# Parents are rows of geno: recurrent %s, donor %s",
            s$rp_name, s$dp_name),
    sprintf("rp_row <- %s", repro_value(s$rp)),
    sprintf("dp_row <- %s", repro_value(s$dp)),
    sprintf("calls_sep <- %s", repro_value(s$calls_sep)),
    "",
    sprintf("# Columns 1 to %d describe the samples; markers start at column %d", ms - 1, ms),
    sprintf("meta <- geno[, colnames(geno)[1:%d]]", ms - 1),
    sprintf("calls <- geno[, colnames(geno)[%d:ncol(geno)]]", ms)
  )

  # Marker QC filters, in the order the app runs them
  if (any_filter) code <- c(code, "", "# Marker QC")
  if (isTRUE(s$apply_par_poly)) {
    code <- c(code, "# Drop markers that are monomorphic between the parents",
              repro_call("parent_poly", c(qc_args, sep_arg), assign = "calls"))
  }
  if (isTRUE(s$apply_par_miss)) {
    code <- c(code, "# Drop markers with a missing parent call",
              repro_suffix(repro_call("parent_missing", qc_args, assign = "calls"), "$par_present"))
  }
  if (isTRUE(s$apply_geno_good)) {
    code <- c(code, "# Drop markers with genotype call errors",
              repro_suffix(repro_call("geno_error",
                                c(qc_args, sep_arg, list(data_type = s$data_type)),
                                assign = "calls"), "$geno_good"))
  }
  if (isTRUE(s$apply_par_homo)) {
    code <- c(code, "# Drop markers where a parent is heterozygous",
              repro_suffix(repro_call("parent_het",
                                c(qc_args, sep_arg, list(na_code = NA)),
                                assign = "calls"), "$par_hom"))
  }
  code <- c(code, if (any_filter) "geno_qc <- cbind(meta, calls)" else "geno_qc <- geno")

  # Map file
  cols <- s$map_cols
  if (identical(s$feedback, "no")) {
    code <- c(
      code, "",
      "# Map file built from the marker names (e.g. S1_778962: chromosome 1, position 778962)",
      repro_call("parse_marker_ns",
                 list(x = repro_sym(sprintf("colnames(geno_qc)[%d:ncol(geno_qc)]", ms)),
                      sep = s$marker_sep, prefix = s$prefix),
                 assign = "map"),
      "map <- order_markers(map)",
      sprintf("map <- map[map[[%s]] %%in%% colnames(calls), ]", repro_value(cols[["snpid"]]))
    )
  } else {
    map_name <- repro_file_name(s$map_name)
    code <- c(
      code, "",
      "# Map file, limited to the markers in the genotype file",
      repro_read_table_code("map", map_name),
      sprintf("common <- intersect(map[[%s]], colnames(geno_qc))", repro_value(s$map_snp_col)),
      sprintf("map <- map[match(common, map[[%s]]), ]", repro_value(s$map_snp_col))
    )
  }

  code <- c(
    code, "",
    "# Genotypes coded as numbers (1 = recurrent, 0 = donor, 0.5 = het)",
    repro_call("proc_kasp",
               list(x = repro_sym("geno_qc"), sample_id = s$sample_id,
                    map_snp_id = cols[["snpid"]], marker_start = ms,
                    chr = cols[["chr"]], chr_pos = cols[["pos"]],
                    kasp_map = repro_sym("map")),
               assign = "geno_num"),
    repro_call("kasp_numeric",
               list(x = repro_sym("geno_num"), rp_row = repro_sym("rp_row"),
                    dp_row = repro_sym("dp_row"), sep = repro_sym("calls_sep"),
                    data_type = s$data_type),
               assign = "geno_num"),
    repro_save_table("cbind(sample = rownames(geno_num), geno_num)", "genotypes_numeric.csv"),
    repro_save_table("map", "map_file.csv")
  )

  files <- stats::setNames(s$data_path, data_name)
  if (!identical(s$feedback, "no")) {
    files <- c(files, stats::setNames(s$map_path, repro_file_name(s$map_name)))
  }

  settings <- s[c("cluster_by", "cluster_focus", "marker_start", "sample_id",
                  "calls_sep", "data_type", "rp", "dp", "rp_name", "dp_name",
                  "apply_par_poly", "apply_par_miss", "apply_geno_good",
                  "apply_par_homo", "feedback")]
  settings$data_file <- data_name
  if (identical(s$feedback, "no")) {
    settings$marker_sep <- s$marker_sep
    settings$prefix <- s$prefix
  } else {
    settings$map_file <- repro_file_name(s$map_name)
  }

  repro_step(
    title = "Read the genotypes, run marker QC and code them as numbers",
    code = code,
    settings = settings,
    files = files
  )
}


#' Step: recurrent parent percentage (RPP) of each progeny
#'
#' @param s List with `rpp` (arguments of calc_rpp_bc() other than `x` and
#'   `map_file`), `plot` (arguments of rpp_barplot() other than `rpp_df`
#'   and `pdf`, or NULL if no plot was drawn), `threshold` (number),
#'   `threshold_note`, `parents` and `replicates`.
#' @noRd
repro_ds_rpp_step <- function(s, width = 8, height = 6) {
  sid <- s$sample_col
  rcol <- s$rpp_col

  code <- c(
    repro_call("calc_rpp_bc",
               c(list(x = repro_sym("geno_num"), map_file = repro_sym("map")), s$rpp),
               assign = "rpp"),
    "",
    paste("# Progenies that reach the RPP threshold", s$threshold_note),
    sprintf("threshold <- %s", repro_value(s$threshold)),
    sprintf("parents <- %s", repro_value(s$parents)),
    sprintf("replicates <- %s  # parent replicates, kept as controls",
            repro_value(as.character(s$replicates))),
    sprintf("rpp_table <- rpp[!rpp[[%s]] %%in%% parents, , drop = FALSE]", repro_value(sid)),
    sprintf("rpp_table$status <- ifelse(rpp_table[[%s]] %%in%% replicates, \"Control\",", repro_value(sid)),
    sprintf("                           ifelse(rpp_table[[%s]] >= threshold, \"Above\", \"Below\"))", repro_value(rcol)),
    sprintf("rpp_table <- rpp_table[order(rpp_table$status == \"Control\", -rpp_table[[%s]]), ]", repro_value(rcol)),
    repro_save_table("rpp_table", "rpp_table.csv")
  )
  if (!is.null(s$plot)) {
    code <- c(
      code, "",
      repro_call("rpp_barplot",
                 c(list(rpp_df = repro_sym("rpp")), s$plot, list(pdf = FALSE)),
                 assign = "rpp_plot"),
      repro_save_plots("list(rpp_plot)", "rpp_barplot.pdf", width, height)
    )
  }

  repro_step(
    title = "Recurrent parent percentage (RPP) of each progeny",
    code = code,
    settings = c(s$rpp, list(threshold = s$threshold, replicates = s$replicates),
                 if (!is.null(s$plot)) list(plot = s$plot)),
    needs = "process"
  )
}


#' Step: heatmap of the progenies
#'
#' @param s List with `fn` (`"cross_qc_heatmap"`, `"cross_qc_heatmap2"` or
#'   `"cross_qc_annotate"`), `args` (its arguments other than `x`,
#'   `map_file` and `pdf`) and, for `cross_qc_annotate`, `zoom_idx` (the
#'   marker columns in the QTL regions).
#' @noRd
repro_ds_heatmap_step <- function(s, width = 10, height = 8) {
  if (identical(s$fn, "cross_qc_annotate")) {
    data_code <- c(
      "# Markers in or nearest to the chosen QTL regions",
      sprintf("zoom_idx <- %s", repro_value(as.integer(s$zoom_idx))),
      "geno_zoom <- geno_num[, zoom_idx, drop = FALSE]",
      "map_zoom <- map[zoom_idx, , drop = FALSE]"
    )
    data_args <- list(x = repro_sym("geno_zoom"), map_file = repro_sym("map_zoom"))
  } else {
    data_code <- character()
    data_args <- list(x = repro_sym("geno_num"), map_file = repro_sym("map"))
  }

  repro_step(
    title = "Heatmap of the progenies",
    code = c(
      data_code,
      repro_call(s$fn, c(data_args, s$args, list(pdf = FALSE)), assign = "heatmaps"),
      repro_save_plots("heatmaps", "heatmaps.pdf", width, height)
    ),
    settings = c(list(plot = s$fn), s$args),
    needs = "process"
  )
}


#' Step: introgression check at the QTLs
#' @param s Arguments of introgression_check() other than `x` and
#'   `map_file`.
#' @noRd
repro_ds_check_step <- function(s) {
  repro_step(
    title = "Introgression check at the QTLs",
    code = c(
      repro_call("introgression_check",
                 c(list(x = repro_sym("geno_num"), map_file = repro_sym("map")), s),
                 assign = "check"),
      repro_save_table("check", "introgression_check.csv")
    ),
    settings = s,
    needs = "process"
  )
}


#' Step: foreground selection
#'
#' @param s List with `geno_name`, `geno_path`, `info_name`, `info_path`,
#'   `sample_col`, `select` (arguments of foreground_select() other than
#'   the data), `upset` (UpSet plot settings, or NULL), `find` (list of
#'   `present` and `absent`, or NULL).
#' @noRd
repro_ds_foreground_step <- function(s, width = 8, height = 6) {
  geno_name <- repro_file_name(s$geno_name)
  info_name <- repro_file_name(s$info_name)
  mcol <- repro_value(s$select$fore_marker_col)

  code <- c(
    repro_read_table_code("geno", geno_name),
    repro_read_table_code("marker_info", info_name),
    "",
    "# Keep the markers that are columns of the genotype file",
    sprintf("marker_info_used <- marker_info[marker_info[[%s]] %%in%% colnames(geno), ]", mcol),
    "",
    "# Line names become row names (missing or repeated names made unique)",
    sprintf("lines <- as.character(geno[[%s]])", repro_value(s$sample_col)),
    "if (anyDuplicated(lines) || anyNA(lines)) {",
    "  lines <- make.unique(ifelse(is.na(lines), \"Unknown\", lines))",
    "}",
    "rownames(geno) <- lines",
    "",
    repro_call("foreground_select",
               c(list(geno_data = repro_sym("geno"),
                      fore_marker_info = repro_sym("marker_info_used")),
                 s$select),
               assign = "fg"),
    "write.csv(data.frame(line = rownames(fg), loci_carried = rowSums(fg == 1, na.rm = TRUE),",
    "                     fg, check.names = FALSE),",
    "          file.path(\"results\", \"foreground_selection.csv\"), row.names = FALSE)"
  )

  if (!is.null(s$upset)) {
    u <- s$upset
    code <- c(
      code, "",
      "# UpSet plot of the favourable-allele combinations",
      sprintf("set_info <- data.frame(sets = marker_info[[%s]], locus = marker_info[[%s]])",
              repro_value(u$marker_col), repro_value(u$locus_col)),
      "upset_plot <- UpSetR::upset(",
      "  fg,",
      "  nsets = ncol(fg),",
      sprintf("  mainbar.y.label = %s,", repro_value(u$mainbar_y_label)),
      sprintf("  sets.x.label = %s,", repro_value(u$sets_x_label)),
      sprintf("  text.scale = %s,", repro_value(u$text_scale)),
      "  set.metadata = list(",
      "    data = set_info,",
      "    plots = list(list(",
      sprintf("      type = \"text\", assign = %s, column = \"locus\",", repro_value(u$assign)),
      sprintf("      colors = rep(%s, length.out = ncol(fg))", repro_value(u$colors)),
      "    ))",
      "  )",
      ")",
      sprintf("pdf(file.path(\"results\", \"upset_plot.pdf\"), width = %s, height = %s, onefile = FALSE)",
              repro_value(width), repro_value(height)),
      "print(upset_plot)",
      "invisible(dev.off())"
    )
  }

  if (!is.null(s$find)) {
    code <- c(
      code, "",
      "# Lines with this combination of favourable alleles",
      repro_call("find_lines",
                 list(mat = repro_sym("fg"),
                      present = as.character(s$find$present),
                      absent = as.character(s$find$absent)),
                 assign = "found_lines"),
      "writeLines(as.character(found_lines), file.path(\"results\", \"found_lines.txt\"))"
    )
  }

  repro_step(
    title = "Foreground selection",
    code = code,
    settings = c(list(genotype_file = geno_name, marker_info_file = info_name,
                      line_column = s$sample_col),
                 s$select,
                 if (!is.null(s$upset)) list(upset = s$upset),
                 if (!is.null(s$find)) list(find_lines = s$find)),
    files = stats::setNames(c(s$geno_path, s$info_path), c(geno_name, info_name))
  )
}
