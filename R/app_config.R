#' Access files in the current app
#'
#' NOTE: If you manually change your package name in the DESCRIPTION,
#' don't forget to update it here and in the config file.
#' For a safer renaming mechanism, use `golem::set_golem_name()`.
#'
#' @param ... character vectors specifying subdirectory and file(s)
#' within your package. Default (none) returns the root of the app.
#'
#' @noRd
app_sys <- function(...) {
  options(shiny.maxRequestSize = 500 * 1024^2)  # Extend file upload size to 500 MB
  system.file(..., package = "panGenomeBreedr")
}




#' Read App Config
#'
#' @param value Value to retrieve from the config file.
#' @param config GOLEM_CONFIG_ACTIVE value. If unset, falls back to R_CONFIG_ACTIVE.
#' If unset, defaults to "default".
#' @param use_parent Logical. If TRUE, scan parent directory for config file.
#' @param file Path to config file. Defaults to internal "golem-config.yml".
#'
#' @noRd
get_golem_config <- function(
    value,
    config = Sys.getenv(
      "GOLEM_CONFIG_ACTIVE",
      Sys.getenv("R_CONFIG_ACTIVE", "default")
    ),
    use_parent = TRUE,
    file = app_sys("golem-config.yml")
) {
  config::get(
    value = value,
    config = config,
    file = file,
    use_parent = use_parent
  )
}


