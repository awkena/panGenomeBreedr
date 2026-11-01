#' Number of Samples per Plate UI Function
#'
#' @description A Shiny module for displaying the number of samples per plate in KASP genotyping data
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel selectInput actionButton icon div
#' @importFrom bslib accordion accordion_panel
#'
mod_mv_nsamples_plate_ui <- function(id) {
  ns <- NS(id)
  ui_sidebar_layout(

    # SIDEBAR ---------------------------------------------------------------

    sidebar = tagList(
      selectizeInput(
        inputId = ns("subset_id"),
        label = label_help(
          "Plate Subset Column",
          "Column that identifies each plate-by-marker subset (the 'plates' column added on upload)."
        ),
        choices = NULL,
        multiple = FALSE
      ),
      selectizeInput(
        inputId = ns("snps_id"),
        label = "SNP ID Column",
        choices = NULL,
        multiple = FALSE
      ),
      selectizeInput(
        inputId = ns("plates_id"),
        label = label_help(
          "Master Plate ID Column",
          "Column with the master plate name. Samples are counted per master plate."
        ),
        choices = NULL,
        multiple = FALSE
      )
    ),

    # MAIN PANEL ------------------------------------------------------------

    tags$p(class = "fw-semibold mb-2", "Samples per plate"),
    DT::DTOutput(ns("nsumm_output"))
  )
}


#' Number of Samples per Plate Server Function
#'
#' @description Server logic for calculating and displaying the number of samples
#'   per plate in KASP genotyping data
#'
#' @param id An ID string that corresponds with the ID used to call the module's UI function
#' @param kasp_data A reactive expression that returns the KASP data frame
#'
#' @return A summary of samples per plate
#'
#' @noRd
#'
#' @importFrom shiny moduleServer reactive reactiveVal observeEvent req updateSelectInput invalidateLater
#'
mod_mv_nsamples_plate_server <- function(id, kasp_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # COLUMN DROPDOWNS ------------------------------------------------------

    # Get colnames from Kasp data
    subset_names <- reactive({
      req(kasp_data())
      return(colnames(kasp_data()))
    })

    # Fill the dropdowns and guess the right column for each
    observe({
      req(subset_names())
      updateSelectizeInput(session,
        server = TRUE,
        inputId = "subset_id",
        choices = subset_names(),
        selected = safe_grep_match(pattern = 'plates',
                                   choices = subset_names())
      )

      updateSelectizeInput(session,
        server = TRUE,
        inputId = "snps_id",
        choices = subset_names(),
        selected = safe_grep_match(pattern = 'snp',
                                   choices = subset_names())
      )

      updateSelectizeInput(session,
        server = TRUE,
        inputId = "plates_id",
        choices = subset_names(),
        selected = safe_grep_match(pattern = 'masterplate',
                                   choices = subset_names())
      )
    })

    # PLATE SUMMARY ---------------------------------------------------------

    # Number of samples per plate
    nsample_plate_result <- reactive({
      req(kasp_data(), input$subset_id, input$snps_id, input$plates_id)
      req(all(c(input$subset_id, input$snps_id, input$plates_id) %in% colnames(kasp_data())))

      tryCatch(
        nsamples_plate(
          x = kasp_data(),
          subset = input$subset_id,
          snp_id = input$snps_id,
          plate_id = input$plates_id
        )$summ,
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Couldn't summarise the plates",
            text = e$message,
            type = "error",
            showCloseButton = TRUE
          )
          NULL
        }
      )
    })

    # Summary table
    output$nsumm_output <- DT::renderDT({
      validate(need(kasp_data(), msg_no_data))
      req(nsample_plate_result())
      df <- data.frame(nsample_plate_result(), check.names = FALSE)
      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
    })

    # Toast once each time the summary changes
    last_summary <- NULL
    observeEvent(nsample_plate_result(), {
      if (identical(nsample_plate_result(), last_summary)) return()
      last_summary <<- nsample_plate_result()
      shinyWidgets::show_toast(
        title = "",
        type = "success",
        text = "Plate summary ready",
        timer = 2000,
        position = "bottom-end"
      )
    })
  })
}

## To be copied in the UI
# mod_mv_nsamples_plate_ui("mv_nsamples_plate_1")

## To be copied in the server
# mod_mv_nsamples_plate_server("mv_nsamples_plate_1")
