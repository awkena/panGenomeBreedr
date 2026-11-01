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
    # An empty numericInput arrives as a logical NA, not a numeric one
    if ((is.numeric(val) || is.logical(val)) && anyNA(val)) {
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


#' Generate a trait input row for genomic coordinates
#'
#' This helper function creates a stylized UI row containing inputs for
#' a trait name, chromosome number, and physical position. It is used
#' within the coordinate definition wizard to allow dynamic addition of multiple loci.
#'
#' @param id A unique integer or character string to suffix the input IDs,
#' ensuring namespace consistency.
#'
#' @return A \code{shiny::tag.list} (div) containing a \code{fluidRow} with
#' \code{textInput} and \code{numericInput} elements.
#'
#' @note The 'Remove' button is programmatically excluded for the first input
#' row (\code{id = 1}) to ensure at least one locus remains defined.
#'
#' @importFrom shiny div fluidRow column textInput numericInput actionButton
#' @noRd
makeTraitInput <- function(id, ns) {
  div(
    id = ns(paste0("trait_", id)),
    style = "border: 1px solid #dee2e6; padding: 15px; margin-bottom: 15px; border-radius: 8px; background-color: #f8f9fa; position: relative;",

    fluidRow(
      # Basic Info
      column(3, textInput(ns(paste0("name_", id)), "Trait/QTL Name:", value = paste0("QTL_", id))),
      column(2, numericInput(ns(paste0("chr_", id)), "Chr:", value = 1, min = 1)),

      # Configuration Checkboxes
      column(7,
             checkboxGroupInput(
               inputId = ns(paste0("type_", id)),
               label = "Annotation Components (Select at least one):",
               choices = c("Point Position" = "pos", "Genomic Range" = "range"),
               selected = "range",
               inline = TRUE
             )
      )
    ),

    fluidRow(
      # Conditional Point Input
      conditionalPanel(
        condition = paste0("input['", ns(paste0("type_", id)), "'].includes('pos')"),
        column(4, numericInput(ns(paste0("pos_", id)), "Exact Position (bp):", value = 1200000))
      ),

      # Conditional Range Inputs
      conditionalPanel(
        condition = paste0("input['", ns(paste0("type_", id)), "'].includes('range')"),
        column(4, numericInput(ns(paste0("start_", id)), "Start (bp):", value = 1000000)),
        column(4, numericInput(ns(paste0("end_", id)), "End (bp):", value = 1400000))
      )
    ),

    # Remove button: only for traits after the first one
    if (id > 1) {
      div(style = "position: absolute; top: 10px; right: 10px;",
          actionButton(ns(paste0("remove_", id)), "", icon = icon("trash"),
                       class = "btn-outline-danger btn-sm", title = "Remove Trait"))
    }
  )
}



#' Render a ghost plate error plot
#'
#' @description This function generates a stylized ggplot2 plot that visually
#' represents a plate layout error. It displays a "ghost" plate background 
#' with an overlayed error message.
#'
#' @param error_msg A character string containing the specific error message to display.
#' @param text_size A numeric value controlling the size of the text elements in the plot.
#' @return A ggplot2 object representing the error plot.
#' @noRd
render_ghost_plate_error <- function(error_msg, text_size = 12) {
  
  # Data setup
  ghost_data <- expand.grid(
    row = factor(rev(LETTERS[1:8])), 
    col = factor(sprintf("%02d", 1:12))
  )

  # Color Palette
  well_grey   <- "#F1F2F6"
  stroke_grey <- "#CED4DA"
  grid_grey   <- "#E9ECEF"
  text_dark   <- "#2D3436"
  error_red   <- "#D63031"

  ggplot2::ggplot(ghost_data, ggplot2::aes(x = col, y = row)) +
    # Draw the empty grid
    ggplot2::geom_tile(fill = "white", color = grid_grey, linewidth = 0.5) +
    
    # Draw the empty wells
    ggplot2::geom_point(
      size = text_size,
      shape = 21,
      fill = well_grey,
      color = stroke_grey,
      stroke = 1,
      alpha = 0.8
    ) +
    
    # Semi-transparent overlay to "de-emphasize" the background
    ggplot2::annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf,
      fill = "white",
      alpha = 0.3
    ) +
    
    # Labels and Theme setup
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(title = "PLATE LAYOUT ERROR") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = text_size,
        color = error_red
      ),
      axis.text = ggplot2::element_text(
        size = text_size,
        face = "bold",
        color = "#B2BEC3"
      )
    ) +
    
    # Central Error Message Label
    ggplot2::annotate(
      "label",
      x = 6.5,
      y = 4.5,
      label = paste0("VALIDATION FAILED:\n", error_msg),
      fill = "white",
      color = text_dark,
      size = text_size * 0.35,
      fontface = "bold",
      label.padding = ggplot2::unit(1, "lines")
    )
}





#' Render a reactable table with custom theming
#'
#' @description This helper function generates a `reactable` table with a default
#'   theme and common display options. It allows for easy and consistent styling
#'   of interactive tables throughout the application.
#'
#' @param data A data frame or tibble to be displayed in the reactable table.
#' @param theme An optional `reactableTheme` object to override the default theme.
#'   If `NULL`, a predefined default theme is used.
#'
#' @return A `reactable` HTML widget displaying the provided data.
#'
#' @noRd
render_reactable <- function(data, theme = NULL) {
  # Add validation check for the data input
  if (is.null(data) || (!is.data.frame(data) && !is.matrix(data))) {
    return(NULL) # Return NULL to avoid error 
  }
  if (is.null(theme)) {
    theme <- reactable::reactableTheme(
      backgroundColor = "hsl(0, 0%, 100%)",
      borderColor = "hsl(0, 0%, 89%)",
      stripedColor = "hsl(0, 0%, 97%)",
      highlightColor = "hsl(0, 0%, 96%)",
      cellPadding = "8px 12px",
      searchInputStyle = list(width = "100%", padding = "8px 12px"),
      headerStyle = list(borderWidth = "1px", padding = "8px 12px")
    )
  }
  reactable::reactable(
    data,
    defaultColDef = reactable::colDef(
      headerClass = "header",
      align = "center",
      minWidth = 100,
      width = 250
    ),
    pagination = TRUE,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(10, 25, 50, 100),
    striped = TRUE,
    highlight = TRUE,
    bordered = TRUE,
    compact = TRUE,
    wrap = TRUE,
    theme = theme,
    style = list(maxWidth = "100%")
  )
}






#' Update sidebar navigation buttons
#'
#' This function dynamically updates the CSS classes of sidebar navigation buttons
#' to visually indicate the active tab. It uses `shinyjs` to toggle Bootstrap 5
#' button classes between outlined and filled states.
#'
#' @param active_id Character string containing the input ID of the currently active button.
#'
#' @return No return value. Called for its side effects on the UI.
#'
#' @noRd
update_sidebar_buttons <- function(active_id) {
  # Updated to match the exact IDs in active_dashboard_view
  buttons <- c(
    "get_db_info",
    "identify_variants_btn",
    "design_kasp_sidebar_btn",
    "pcil_selection_btn"
  )

  for (btn_id in buttons) {
    if (btn_id == active_id) {
      shinyjs::removeClass(btn_id, "btn-outline-primary")
      shinyjs::addClass(btn_id, "btn-primary")
    } else {
      shinyjs::removeClass(btn_id, "btn-primary")
      shinyjs::addClass(btn_id, "btn-outline-primary")
    }
  }
}


#' Path to an example data file shipped with panGenomeBreedr
#'
#' @param ... Path components under `inst/extdata`, e.g. `"Genotyping_141.010_01.csv"`.
#' @return The full path to the file.
#' @noRd
example_file <- function(...) {
  system.file("extdata", ..., package = "panGenomeBreedr", mustWork = TRUE)
}


#' "Load example data" link placed under a file upload
#'
#' @param id The namespaced input ID of the link (e.g. `ns("load_example")`).
#' @param label Link text.
#' @return A Shiny tag.
#' @noRd
example_data_link <- function(id, label = "Load example data") {
  div(
    style = "margin-top: -10px; margin-bottom: 12px;",
    actionLink(id, label = tagList(icon("flask"), label), class = "small")
  )
}


#' Badge shown while a tab is using example data
#'
#' @param file_name The example file being used.
#' @param clear_id Optional namespaced input ID for a "clear" link.
#' @param note Optional extra text or tags, e.g. a download link for a file
#'   that is not shipped with the package.
#' @return A Shiny tag.
#' @noRd
example_data_badge <- function(file_name, clear_id = NULL, note = NULL) {
  div(
    class = "small",
    style = paste(
      "margin-bottom: 12px; padding: 6px 10px; border-radius: 6px;",
      "background: #eef6fb; border: 1px solid #cfe3f1; color: #1e3a5f;"
    ),
    icon("flask"), " Using example data: ", tags$code(file_name),
    if (!is.null(clear_id)) tagList(" \u00b7 ", actionLink(clear_id, "clear")),
    if (!is.null(note)) div(style = "margin-top: 6px;", note)
  )
}


# ---------------------------------------------------------------------------
# Layout building blocks (used by the KASP Marker Validation modules)
# ---------------------------------------------------------------------------

#' Settings sidebar next to the results
#'
#' A collapsible bslib sidebar for the inputs, with the results beside it.
#'
#' @param sidebar Tags for the sidebar (inputs).
#' @param ... Tags for the main area (results).
#' @param width Sidebar width in pixels.
#' @return A Shiny tag.
#' @noRd
ui_sidebar_layout <- function(sidebar, ..., width = 350) {
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(sidebar, width = width, open = "desktop"),
    border = FALSE,
    fillable = FALSE,
    ...
  )
}


#' Input label with a hover hint
#'
#' @param label Label text.
#' @param help Plain-language explanation shown on hover.
#' @return A Shiny tag to pass as an input's `label`.
#' @noRd
label_help <- function(label, help) {
  tags$span(
    label,
    bslib::tooltip(
      tags$span(
        icon("circle-info"),
        class = "text-muted ms-1",
        style = "cursor: help; font-size: 0.85em;"
      ),
      help,
      placement = "right"
    )
  )
}


#' Results card with a title and optional header buttons
#'
#' @param title Card title.
#' @param ... Card body.
#' @param actions Tags shown on the right of the header (e.g. buttons).
#' @param full_screen Show the expand-to-full-screen button?
#' @return A Shiny tag.
#' @noRd
ui_result_card <- function(title, ..., actions = NULL, full_screen = FALSE) {
  bslib::card(
    full_screen = full_screen,
    bslib::card_header(
      class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
      tags$span(class = "fw-semibold", title),
      if (!is.null(actions)) div(class = "d-flex gap-2", actions)
    ),
    bslib::card_body(...)
  )
}


#' Gear button that opens plot appearance options
#'
#' Do not use `label_help()` inside: a tooltip nested in a popover throws
#' JavaScript errors in bslib 0.6.
#'
#' @param ... Inputs to show in the popover.
#' @param title Popover title.
#' @return A Shiny tag.
#' @noRd
plot_options_button <- function(..., title = "Plot appearance") {
  bslib::popover(
    tags$button(
      type = "button",
      class = "btn btn-sm btn-outline-secondary",
      icon("sliders"), " Appearance"
    ),
    ...,
    title = title,
    placement = "bottom"
  )
}


#' Download button that opens the file options
#'
#' Creates the inputs `<id>_name`, `<id>_type`, `<id>_width` and
#' `<id>_height`; the download handler itself is `output[[id]]`.
#'
#' @param ns The module's namespace function.
#' @param id Output ID of the download handler.
#' @param file_name Default file name, without extension.
#' @param formats Named vector of file formats.
#' @param width,height Default size in inches.
#' @param note Optional text shown above the save button.
#' @return A Shiny tag.
#' @noRd
plot_download_button <- function(ns, id, file_name,
                                 formats = c("PDF" = "pdf", "PNG" = "png"),
                                 width = 8, height = 6, note = NULL) {
  bslib::popover(
    tags$button(
      type = "button",
      class = "btn btn-sm btn-success",
      icon("download"), " Download"
    ),
    textInput(ns(paste0(id, "_name")), "File name", value = file_name),
    selectInput(ns(paste0(id, "_type")), "Format", choices = formats),
    bslib::layout_columns(
      numericInput(ns(paste0(id, "_width")), "Width (in)", value = width, min = 1, max = 30),
      numericInput(ns(paste0(id, "_height")), "Height (in)", value = height, min = 1, max = 30)
    ),
    if (!is.null(note)) tags$p(class = "small text-muted", note),
    downloadButton(ns(id), "Save", class = "btn-success w-100"),
    title = "Download plot",
    placement = "bottom"
  )
}


#' Main action button
#'
#' Shows a spinner and blocks double clicks while the server is busy
#' (bslib >= 0.7.0); falls back to a plain action button otherwise.
#'
#' @param id Input ID.
#' @param label Button text.
#' @param icon_name Font Awesome icon name.
#' @return A Shiny tag.
#' @noRd
run_button <- function(id, label, icon_name = "play") {
  if (utils::packageVersion("bslib") >= "0.7.0") {
    bslib::input_task_button(
      id, label,
      icon = icon(icon_name),
      label_busy = "Working...",
      type = "success",
      class = "w-100 mt-2"
    )
  } else {
    actionButton(id, label, icon = icon(icon_name), class = "btn-success w-100 mt-2")
  }
}


#' Value of an expression, or NULL if it is not ready
#'
#' Turns a `req()` stop or an error into `NULL`, so an output can show a
#' placeholder message with `validate(need(...))` instead of staying blank.
#'
#' @param expr Expression, usually a reactive call.
#' @return The value or NULL.
#' @noRd
if_ready <- function(expr) {
  tryCatch(expr, error = function(e) NULL)
}


#' Placeholder messages shown in empty results
#' @noRd
msg_no_data <- "Upload a KASP file, or load the example data, to see results here."

#' @noRd
msg_no_color <- "Results appear here once the genotype calls are colour-coded (see Color Code Calls)."
