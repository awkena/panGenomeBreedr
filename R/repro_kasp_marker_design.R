#' Reproducible-script step for the KASP Marker Design tab
#'
#' Takes the settings of the last successful design run, as a plain list,
#' and returns a `repro_step()` (see `R/utils_repro.R`).
#'
#' @name repro_marker_design
#' @keywords internal
#' @noRd
NULL


#' Header title of the KASP marker design script
#' @noRd
repro_md_title <- "panGenomeBreedr KASP marker design"

#' Script file name for the KASP Marker Design tab
#' @noRd
repro_md_script <- "kasp_marker_design.R"


#' Step: design KASP markers and make the Intertek table
#'
#' @param s List with:
#'   - `source`: `"vcf"` or `"matrix"`.
#'   - `variant_file`, `variant_name`: local path and original name of the
#'     VCF or genotype matrix.
#'   - `genome_file`, `genome_name`: local path and original name of the
#'     reference genome.
#'   - `genome_url`: optional URL the script downloads the genome from when
#'     it is not in `data/` (used for the example data).
#'   - `design`: arguments of kasp_marker_design() other than the files.
#'   - `intertek`: arguments of make_intertek_table() other than
#'     `marker_data`.
#' @noRd
repro_md_step <- function(s) {
  variant_name <- repro_file_name(s$variant_name)
  genome_name <- repro_file_name(s$genome_name)
  variant_path <- repro_data_path(variant_name)
  genome_path <- repro_data_path(genome_name)

  # Variant data: a VCF path, or a genotype matrix read into gt_df
  if (identical(s$source, "vcf")) {
    read_code <- character()
    files_args <- list(vcf_file = variant_path, gt_df = NULL)
  } else {
    is_excel <- grepl("\\.xlsx?$", variant_name, ignore.case = TRUE)
    read_code <- c(
      "# Genotype matrix",
      if (is_excel) {
        sprintf("gt_df <- as.data.frame(readxl::read_excel(%s, .name_repair = \"minimal\"))",
                repro_value(variant_path))
      } else {
        sprintf("gt_df <- read.csv(%s, stringsAsFactors = FALSE)", repro_value(variant_path))
      },
      ""
    )
    files_args <- list(vcf_file = NULL, gt_df = repro_sym("gt_df"))
  }

  # Reference genome, downloaded if it is not in data/ and a link is known
  genome_code <- c(
    "# Reference genome",
    sprintf("genome_file <- %s", repro_value(genome_path)),
    if (!is.null(s$genome_url)) {
      c("if (!file.exists(genome_file)) {",
        sprintf("  download.file(%s, genome_file, mode = \"wb\")", repro_value(s$genome_url)),
        "}")
    },
    ""
  )

  design_args <- c(
    files_args,
    s$design,
    list(
      genome_file = repro_sym("genome_file"),
      # Writes one alignment_<marker>.pdf per marker
      save_alignment = TRUE,
      plot_file = "results"
    )
  )

  repro_step(
    title = "Design KASP markers and make the Intertek table",
    code = c(
      read_code,
      genome_code,
      repro_call("kasp_marker_design", design_args, assign = "markers"),
      "",
      "# Assay table to send to Intertek",
      repro_call("make_intertek_table",
                 c(list(marker_data = repro_sym("markers")), s$intertek),
                 assign = "intertek_table"),
      "write.csv(intertek_table, file.path(\"results\", \"KASP_markers_for_Intertek.csv\"), row.names = FALSE)"
    ),
    settings = c(
      list(variant_data = s$source, variant_file = variant_name,
           genome_file = genome_name),
      s$design, s$intertek
    ),
    files = stats::setNames(c(s$variant_file, s$genome_file),
                            c(variant_name, genome_name))
  )
}
