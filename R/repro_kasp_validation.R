#' Reproducible-script steps for the KASP marker validation modules
#'
#' Each function takes the settings a module used, as a plain list, and
#' returns a `repro_step()`. The modules call them from a reactive and
#' register the result with a `repro_recorder()` (see `R/utils_repro.R`).
#'
#' Objects created by the script, for the steps that build on each other:
#' `kasp` (read_kasp_csv() result), `kasp_data` (the Data section with a
#' `plates` column) and `kasp_coded` (kasp_color() result).
#'
#' @name repro_kasp
#' @keywords internal
#' @noRd
NULL


#' Header title of the KASP validation scripts
#' @noRd
repro_kasp_title <- "panGenomeBreedr KASP marker validation"


#' Script file name for the KASP Marker Validation tab
#' @noRd
repro_kasp_script <- "marker_validation.R"


#' Plain names of the KASP steps, shown when a step still has to run
#' @noRd
repro_kasp_labels <- c(
  read = "Upload Data: choose a KASP file and click Run Analysis",
  call_rate = "Data Summary: marker call rates",
  nsamples = "Data Summary: plate statistics",
  color = "Color Code Calls: pick the subset and genotype call columns",
  qc = "QC Plots: pick the FAM, HEX and call columns",
  pred = "Prediction Summary: pick the SNP, group and call columns"
)


#' Code lines that save a list of plots to one multi-page PDF
#' @noRd
repro_save_plots <- function(obj, file, width = 8, height = 6) {
  c(
    sprintf("pdf(file.path(\"results\", %s), width = %s, height = %s, onefile = TRUE)",
            repro_value(file), repro_value(width), repro_value(height)),
    sprintf("for (p in %s) print(p)", obj),
    "invisible(dev.off())"
  )
}


#' Code line that saves a table to results/
#' @noRd
repro_save_table <- function(expr, file) {
  sprintf("write.csv(%s, file.path(\"results\", %s), row.names = FALSE)",
          expr, repro_value(file))
}


#' Step: read the KASP file and add the plates column
#'
#' @param s List with `name` (original file name), `path` (local copy),
#'   `datatype`, `row_tags` (character vector) and `spacing`.
#' @noRd
repro_kasp_read_step <- function(s) {
  name <- repro_file_name(if (is.null(s$name)) "kasp_results.csv" else s$name)

  args <- list(file = repro_data_path(name))
  if (identical(s$datatype, "raw")) {
    args$row_tags <- s$row_tags
    args$spacing <- s$spacing
  }
  args$data_type <- s$datatype

  repro_step(
    title = "Read the KASP results file",
    code = c(
      repro_call("read_kasp_csv", args, assign = "kasp"),
      "",
      "# Genotype calls, with a 'plates' column (MasterPlate_SNPID) for each",
      "# plate-by-marker subset",
      "kasp_data <- as.data.frame(kasp$Data)",
      "if (!\"plates\" %in% colnames(kasp_data)) {",
      "  kasp_data$plates <- paste(kasp_data$MasterPlate, kasp_data$SNPID, sep = \"_\")",
      "}",
      repro_save_table("kasp_data", "kasp_data.csv")
    ),
    settings = list(
      file = name,
      data_type = s$datatype,
      row_tags = if (identical(s$datatype, "raw")) s$row_tags,
      spacing = if (identical(s$datatype, "raw")) s$spacing
    ),
    files = stats::setNames(s$path, name)
  )
}


#' Step: call rates per marker and plate
#' @param s Arguments of kasp_call_rate() other than `x`.
#' @noRd
repro_kasp_call_rate_step <- function(s, width = 8, height = 6) {
  repro_step(
    title = "Call rates per marker and plate",
    code = c(
      repro_call("kasp_call_rate", c(list(x = repro_sym("kasp_data")), s),
                 assign = "call_rates"),
      repro_save_table("call_rates$marker_summary", "call_rate_by_marker.csv"),
      repro_save_table("call_rates$plate_summary", "call_rate_by_plate.csv"),
      repro_save_plots(
        "call_rates[c(\"marker_call_rate\", \"no_call_reason\", \"plate_call_rate\")]",
        "call_rate_plots.pdf", width, height
      )
    ),
    settings = s,
    needs = "read"
  )
}


#' Step: number of samples per plate
#' @param s Arguments of nsamples_plate() other than `x`.
#' @noRd
repro_kasp_nsamples_step <- function(s) {
  repro_step(
    title = "Samples per plate",
    code = c(
      repro_call("nsamples_plate", c(list(x = repro_sym("kasp_data")), s),
                 assign = "plate_counts"),
      repro_save_table("plate_counts$summ", "samples_per_plate.csv")
    ),
    settings = s,
    needs = "read"
  )
}


#' Step: colour-code the genotype calls
#' @param s Arguments of kasp_color() other than `x`.
#' @noRd
repro_kasp_color_step <- function(s) {
  repro_step(
    title = "Colour-code the genotype calls per plate",
    code = c(
      repro_call("kasp_color", c(list(x = repro_sym("kasp_data")), s),
                 assign = "kasp_coded")
    ),
    settings = s,
    needs = "read"
  )
}


#' Step: QC cluster plots and plate layouts
#' @param s Arguments of kasp_qc_ggplot2() other than `x` and `pdf`.
#' @param layout Arguments of plot_plate() other than `x`, `filename` and
#'   `pdf`: `well`, `geno_call` and `snp_id`.
#' @noRd
repro_kasp_qc_step <- function(s, layout, width = 8, height = 6) {
  repro_step(
    title = "QC cluster plots and plate layouts",
    code = c(
      repro_call("kasp_qc_ggplot2",
                 c(list(x = repro_sym("kasp_coded"), pdf = FALSE), s),
                 assign = "qc_plots"),
      repro_save_plots("qc_plots", "qc_plots.pdf", width, height),
      "",
      repro_call("plot_plate",
                 c(list(x = repro_sym("kasp_coded")), layout,
                   list(filename = NULL, pdf = FALSE)),
                 assign = "plate_layouts"),
      repro_save_plots("plate_layouts", "plate_layouts.pdf", width, height)
    ),
    settings = c(s, list(plate_layout = layout)),
    needs = "color"
  )
}


#' Step: prediction summary of the positive controls
#' @param s Arguments of pred_summary() other than `x`.
#' @param plot Arguments of pred_summary_plot() other than `x`, or NULL to
#'   leave out the plots.
#' @noRd
repro_kasp_pred_step <- function(s, plot = NULL, width = 8, height = 6) {
  code <- c(
    repro_call("pred_summary", c(list(x = repro_sym("kasp_coded")), s),
               assign = "pred"),
    repro_save_table("pred$summ", "prediction_summary.csv"),
    "for (plate in names(pred$plates)) {",
    "  write.csv(pred$plates[[plate]],",
    "            file.path(\"results\", paste0(\"prediction_status_\", plate, \".csv\")),",
    "            row.names = FALSE)",
    "}"
  )
  if (!is.null(plot)) {
    code <- c(
      code, "",
      repro_call("pred_summary_plot",
                 c(list(x = repro_sym("as.data.frame(pred$summ)")), plot),
                 assign = "pred_plots"),
      repro_save_plots("pred_plots", "prediction_plots.pdf", width, height)
    )
  }
  repro_step(
    title = "Prediction summary of the positive controls",
    code = code,
    settings = c(s, if (!is.null(plot)) list(plot = plot)),
    needs = "color"
  )
}
