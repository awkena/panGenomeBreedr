#' Check a set of Shiny inputs before running an action, and show a modal
#' naming whichever ones are still empty.
#'
#' Meant to replace a bare `req(...)` guard at the top of a button's
#' `observeEvent()`, where a missing input currently just cancels the handler
#' silently. Instead, this tells the user exactly what they still need to
#' fill in.
#'
#' @param fields A named list; each name is the human-readable label shown to
#'   the user, and each value is the current value of that field (e.g.
#'   `input$some_field`, or a piece of reactive state the action depends on).
#' @param title Modal title.
#'
#' @return `TRUE` (invisibly) if every field is filled in. If any are empty,
#'   shows a modal listing them and returns `FALSE`, so the caller can
#'   `if (!validate_required_inputs(...)) return()`.
#'
#' @noRd
validate_required_inputs <- function(fields, title = "Missing Information") {
  is_missing <- function(val) {
    if (is.null(val)) {
      return(TRUE)
    }
    if (is.data.frame(val)) {
      return(nrow(val) == 0)
    }
    if (length(val) == 0) {
      return(TRUE)
    }
    if (is.character(val) && all(trimws(val) == "")) {
      return(TRUE)
    }
    if (is.numeric(val) && anyNA(val)) {
      return(TRUE)
    }
    FALSE
  }

  missing_fields <- names(fields)[vapply(fields, is_missing, logical(1))]

  if (length(missing_fields) == 0) {
    return(invisible(TRUE))
  }

  shiny::showModal(shiny::modalDialog(
    title = shiny::div(
      class = "d-flex align-items-center",
      shiny::icon(
        "triangle-exclamation",
        class = "text-warning me-2",
        style = "font-size: 1.5rem;"
      ),
      shiny::tags$b(title)
    ),
    size = "m",
    easyClose = TRUE,
    fade = TRUE,
    shiny::div(
      class = "text-muted mb-2",
      "Please fill in the following before continuing:"
    ),
    shiny::tags$ul(lapply(missing_fields, shiny::tags$li)),
    footer = shiny::modalButton("Close")
  ))

  FALSE
}
