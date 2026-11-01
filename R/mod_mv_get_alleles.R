#' MV Get Alleles UI Function
#'
#' @description A Shiny module for getting alleles.
#'
#' @param id Internal parameter for {shiny}.
#'
#'
#' @noRd
#'
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel selectInput textInput tableOutput
#' @importFrom bslib accordion accordion_panel
#'
mod_mv_get_alleles_ui <- function(id) {
  ns <- NS(id)
  ui_sidebar_layout(

    # SIDEBAR ---------------------------------------------------------------

    sidebar = tagList(
      selectInput(
        inputId = ns("plates_pres"),
        label = "Plate",
        choices = NULL,
        multiple = FALSE
      ),
      selectInput(
        inputId = ns("data_type_id"),
        label = "Genotyping Platform",
        choices = c("kasp", "agriplex"),
        multiple = FALSE
      ),
      textInput(
        inputId = ns("sep_id"),
        label = label_help(
          "Genotype Separator",
          "Character between the two alleles in a call, e.g. ':' in A:G."
        ),
        value = ":"
      )
    ),

    # MAIN PANEL ------------------------------------------------------------

    bslib::layout_columns(
      col_widths = c(5, 7),
      # Alleles table
      div(
        tags$p(class = "fw-semibold mb-2", "Identified alleles"),
        DT::DTOutput(ns("alleles_output"))
      ),
      # Genotypes table
      div(
        tags$p(class = "fw-semibold mb-2", "Genotype classes"),
        DT::DTOutput(ns("genotypes_output"))
      )
    )
  )
}


#' MV Get Alleles Server Functions
#'
#' @description Server logic for the get alleles module
#'
#' @param id Internal parameter for {shiny}.
#' @param kasp_data A reactive value  returned from read_kasp_csv
#' @param input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shiny moduleServer NS req reactive observeEvent
#'
#'
#' @noRd
mod_mv_get_alleles_server <- function(id, kasp_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # PLATE DROPDOWN --------------------------------------------------------

    # Plates in the data
    distinct_plates <- reactive({
      req(kasp_data())
      uniq_plates(kasp_data())
    })

    # Fill the plate dropdown
    observe({
      req(distinct_plates())
      updateSelectInput(session,
        inputId = "plates_pres",
        choices = distinct_plates()
      )
    })

    # ALLELES ---------------------------------------------------------------

    # Calls for the selected plate
    Get_calls <- reactive({
      req(kasp_data(), input$plates_pres)
      get_calls(x = kasp_data(), a = input$plates_pres)
    })

    # Alleles and genotypes from the calls
    obtain_alleles <- reactive({
      req(Get_calls(), input$data_type_id, input$sep_id)
      get_alleles(
        x = Get_calls(),
        data_type = input$data_type_id,
        sep = input$sep_id
      )
    })

    # OUTPUTS ---------------------------------------------------------------

    # Alleles table
    output$alleles_output <- DT::renderDT({
      validate(need(kasp_data(), msg_no_data))
      req(obtain_alleles())
      alleles_df(x = obtain_alleles()$alleles)
    })

    # Genotypes table
    output$genotypes_output <- DT::renderDT({
      validate(need(kasp_data(), "Nothing to show yet."))
      req(obtain_alleles())
      genotypes(obtain_alleles()$genotypes)
    })

    # Toast when a new plate is shown
    observeEvent(obtain_alleles(), {
      shinyWidgets::show_toast(
        title = '',
        text = paste("Showing alleles for", input$plates_pres),
        type = "success",
        timer = 2000,
        position = "bottom-end", width = "30%"
      )
    })
  })
}

## To be copied in the UI
# mod_mv_get_alleles_ui("mv_get_alleles_1")

## To be copied in the server
# mod_mv_get_alleles_server("mv_get_alleles_1")
