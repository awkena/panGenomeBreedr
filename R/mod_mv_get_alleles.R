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
#' @importFrom DT DTOutput
#'
mod_mv_get_alleles_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Get alleles
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = ns("data_type_id"),
          label = "Select Data Source",
          choices = c("kasp", "agriplex"),
          multiple = FALSE
        ),
        textInput(
          inputId = ns("sep_id"),
          label = "Genotype Call Separator",
          value = ":"
        )
      ),
      mainPanel(
        # Allele output
        bslib::accordion(
          width = "100%",
          open = "Identified Alleles",
          bslib::accordion_panel(
            "Identified Alleles",
            selectInput(ns("plates_pres"),
                        label = "Choose plate",
                        choices = NULL,
                        multiple = FALSE,
                        width = '40%'
            ),
            DT::DTOutput(ns("alleles_output"))
          ) , style = "margin-bottom: 15px;"
        ),
        # Genotype output
        bslib::accordion(
          width = "100%",
          open = "Genotype Classifications",
          bslib::accordion_panel(
            "Genotype Classifications",
            DT::DTOutput(ns("genotypes_output"))
          )
        )
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
#' @importFrom DT renderDT
#'
#' @noRd
mod_mv_get_alleles_server <- function(id,kasp_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    distinct_plates <- reactive({
      req(kasp_data)
      uniq_plates(kasp_data)
    })

    # Update for get alleles entry for user.
    observe({
      req(distinct_plates())
      updateSelectInput(session,
                        inputId = "plates_pres",
                        choices = distinct_plates()
      )
    })

    # Get genotype calls
    Get_calls <- reactive({
      req(kasp_data, input$plates_pres)
      get_calls(x = kasp_data, a = input$plates_pres)
    })

    # Get alleles
    obtain_alleles <- reactive({

      req(Get_calls(), input$data_type_id , input$sep_id)
      get_alleles(x = Get_calls(),
                                   data_type = input$data_type_id,
                                   sep = input$sep_id)
    })

    # Display alleles from get_alleles()
    observe({
      req(obtain_alleles())

      output$alleles_output <- renderDT({
        req(obtain_alleles())
       alleles_df(x = obtain_alleles()$alleles)
      })

      # Display genotypes from get_alleles()
      output$genotypes_output <- renderDT({
        req(obtain_alleles())
        genotypes(obtain_alleles()$genotypes)
      })

      # show toast
      shinyWidgets::show_toast(
        title = 'Showing Alleles for',
        text = input$plates_pres,
        type = 'success',
        timer = 2000,
        position = "bottom-end",width = '30%'
      )

    })

  })
}

## To be copied in the UI
# mod_mv_get_alleles_ui("mv_get_alleles_1")

## To be copied in the server
# mod_mv_get_alleles_server("mv_get_alleles_1")
