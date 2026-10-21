#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' @param auto_install Logical. If `TRUE` (default) and the session is
#'   interactive, offers to install missing packages the app needs, after
#'   asking for confirmation. Nothing is installed without consent. If `FALSE`,
#'   or in a non-interactive session, missing packages stop the app with
#'   install instructions instead.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
    onStart = NULL,
    auto_install = TRUE,
    options = list(launch.browser = TRUE),
    enableBookmarking = NULL,
    uiPattern = "/",
    ...) {
  # CRAN packages (listed in Suggests) that the app UI and server call directly
  app_packages <- c(
    "colourpicker", "config", "data.table", "DT", "future", "leaflet",
    "openxlsx", "promises", "reactable", "readxl", "shinyalert", "shinybusy",
    "shinycssloaders", "shinyFiles", "shinyjs", "shinyWidgets", "UpSetR",
    "vcfR", "writexl"
  )

  # Bioconductor packages needed only by the KASP Marker Design tab
  bioc_packages <- c("Biostrings", "BSgenome", "GenomicRanges", "IRanges", "msa")

  is_missing <- function(pkgs) {
    pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  }

  missing_app <- is_missing(app_packages)

  if (length(missing_app) > 0) {
    install_cmd <- sprintf("install.packages(c(%s))",
                           paste(sprintf("'%s'", missing_app), collapse = ", "))

    consent <- auto_install && interactive() &&
      isTRUE(utils::askYesNo(
        sprintf("The panGenomeBreedr app needs these packages: %s\nInstall them now?",
                paste(missing_app, collapse = ", "))
      ))

    if (consent) {
      utils::install.packages(missing_app)
      missing_app <- is_missing(missing_app)
    }

    if (length(missing_app) > 0) {
      stop(
        sprintf("Missing packages required by the app: %s\nInstall with: %s",
                paste(missing_app, collapse = ", "), install_cmd),
        call. = FALSE
      )
    }
  }

  # Bioconductor packages do not block launch; only marker design needs them
  missing_bioc <- is_missing(bioc_packages)

  if (length(missing_bioc) > 0) {
    message(
      "The KASP Marker Design tab needs these Bioconductor packages: ",
      paste(missing_bioc, collapse = ", "), "\n",
      "Install with:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
      sprintf("  BiocManager::install(c(%s))",
              paste(sprintf("'%s'", missing_bioc), collapse = ", "))
    )
  }

  # A humour message
   cat('================================================\n')
   cat("panGenomeBreedr is live!\U1F680,Let's breed smarter!\n")
   cat('================================================\n')
   # force open in web browser
   options(shiny.launch.browser = TRUE)


  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
