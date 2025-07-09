#' ds_foreground_select UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_ds_foreground_select_ui <- function(id) {
  ns <- NS(id)
  tagList(

    navset_card_underline(
      nav_panel(
        title = "Trait Introgression Hypothesis Testing",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
          # Side
          ),

          # Main content area
          bslib::input_switch(id = ns("configure"), label = "Configure Plot", value = FALSE),
          conditionalPanel(
            condition = paste0('input["', ns("configure"), '"] == true'),
            fluidRow(
              # Card 1: Heatmap Configuration
              column(
                width = 4,
                bslib::card(
                  max_height = "850px",
                  height = "100%",
                  bslib::card_header("Heatmap Configuration", class = "bg-success"),
                  bslib::card_body(
                    # Upset parameters
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

#' ds_foreground_select Server Functions
#'
#' @noRd
mod_ds_foreground_select_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

  })
}


## To be copied in the UI
# mod_ds_foreground_select_ui("ds_foreground_select_1")

## To be copied in the server
# mod_ds_foreground_select_server("ds_foreground_select_1")

