#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' @param auto_install installs used packages not found in namespace.
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
  # Check all suggested packages at the start
  suggested_packages <- c(
    "DT", "data.table", "fontawesome", "openxlsx", "reactable",
    "readxl", "shinyWidgets", "shinyalert", "shinybusy", "shinyjs",
    "vcfR", "writexl",'UpSetR', 'shinyFiles'
  )

  missing_packages <- suggested_packages[!sapply(suggested_packages,
                                                 requireNamespace,
                                                 quietly = TRUE
  )]

  if (length(missing_packages) > 0) {
    if (auto_install) {
      message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
      install.packages(missing_packages)
    } else {
      stop(
        sprintf("Missing packages: %s\nInstall with: install.packages(c(%s))",
                paste(missing_packages, collapse = ", "),
                paste(sprintf("'%s'", missing_packages), collapse = ", ")),
        call. = FALSE
      )
    }
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
