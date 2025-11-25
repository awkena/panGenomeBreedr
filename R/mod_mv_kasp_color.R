#' KASP Color UI Function
#'
#' @description A shiny Module for Kasp_color
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel selectInput textInput actionButton icon div
#' @importFrom bslib accordion accordion_panel
#'
mod_mv_kasp_color_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = ns("Subset_names"),
          label = "Select Plate Subsetting Column",
          choices = NULL, multiple = FALSE
        ),
        selectInput(
          inputId = ns("geno_call_col"),
          label = "Select Genotype Call Column",
          choices = NULL, multiple = FALSE
        ),
        textInput(
          inputId = ns("sep"),
          label = "Specify Genotype Call Separator",
          value = ":"
        ),
        textInput(
          inputId = ns("uncallable"),
          label = "List Uncallable Genotype Values",
          value = "Uncallable"
        ),
        textInput(
          inputId = ns("unused"),
          label = "List Unused Genotype Values",
          value = "?"
        ),
        textInput(
          inputId = ns("blank"),
          label = "Specify No Template Control Calls",
          value = "NTC"
        ),
        textInput(
          inputId = ns("others"),
          label = "List Non-genotype Call Values",
          value = "Missing , Bad , Dupe , Over , Short"
        ),
        selectInput(
          inputId = ns("color_choose"),
          label = "Colors for FAM, HEX, Het (in order)",
          choices = grDevices::colors(),
          multiple = TRUE,
          selected = c("blue", "gold", "forestgreen")
        )
      ),
      mainPanel(
        bslib::accordion(
          open = TRUE,
          bslib::accordion_panel(
            title = "Plate Status After Applying FAM & HEX Colors",
            DT::DTOutput(outputId = ns("kasp_color_code_stat"))
          ),
          style = "margin-bottom: 15px;"
        ),
        bslib::accordion(
          bslib::accordion_panel(
            title = "Color-Coded KASP Genotype Calls",
            selectInput(
              inputId = ns("subplates"),
              choices = NULL,
              label = "Select Plate", multiple = FALSE
            ),
            DT::DTOutput(outputId = ns("kasp_color_code"))
          ),
          style = "margin-bottom: 15px;"
        )
      )
    )
  )
}
#' KASP Color Server Functions
#'
#' @param id The module ID
#' @param kasp_data The KASP genotyping data from read_kasp_csv
#'
#' @noRd
#'
#' @importFrom shiny moduleServer req observeEvent reactiveVal reactive
#' @importFrom shiny updateSelectInput
#'
mod_mv_kasp_color_server <- function(id, kasp_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # get colnames for drop down
    sub_names <- reactive({
      req(kasp_data)
      colnames(kasp_data)
    })

    # Populate choices in subset
    observe({
      req(sub_names())
      updateSelectInput(session,
                        inputId = "Subset_names",
                        choices = sub_names(),
                        selected = safe_grep_match(pattern = 'plates',
                                                   choices =sub_names())
      )

      updateSelectInput(session,
                        inputId = "geno_call_col",
                        choices = sub_names(),
                        selected = safe_grep_match(pattern = 'call',
                                                   choices = sub_names())
                        )
    })

    # Color code.
    color_coded <- reactive({
      req(
        kasp_data,
        input$Subset_names,
        input$geno_call_col,
        length(input$color_choose) >= 3,  # length of color must be greater than or equal to 3
        input$uncallable,
        input$unused,
        input$blank,
        input$others,
        input$sep
      )

      tryCatch({
        kasp_color(
          x = kasp_data,
          uncallable = trimws(input$uncallable), # remove whitespaces
          unused = trimws(input$unused),# remove whitespaces
          blank = trimws(input$blank), # remove whitespaces
          others = trimws(unlist(strsplit(input$others,split =  ","))),
          sep = check_sep(input$sep), # custom function to check separator
          subset = input$Subset_names,
          geno_call = input$geno_call_col,
          assign_cols = c(
            FAM = input$color_choose[1],
            HEX = input$color_choose[2],
            het = input$color_choose[3]
          )
        )
      }, error = function(e) {
        shinyWidgets::show_alert("Error", paste("Color coding failed:", e$message), type = "error")
        NULL
      })
    })

    # Output: Plate names
    names_col <- reactive({
      req(color_coded())
      names(color_coded())
    })

    # Populate dropdown for subplates
    observe({
      req(names_col())
      updateSelectInput(session, inputId = "subplates", choices = names_col())
    })

    # Display color-coded KASP data
    observe({
      req(input$subplates, color_coded())

      output$kasp_color_code <- DT::renderDT({
        DT::datatable(as.data.frame(color_coded()[[input$subplates]]),
                      options = list(pageLength = 10, scrollX = TRUE)
        )
      })

      # Show toast.
      shinyWidgets::show_toast(
        title = "",
        type = "success",
        text = 'Genotype calls color-coded successfully',
        timer = 2000,
        position = "bottom-end"
      )
    })

    # Status of plates
    color_coded_stat <- reactive({
      req(kasp_data, input$Subset_names, input$geno_call_col)

      kasp_color_stat.give(
        x = kasp_data,
        subset = input$Subset_names,
        geno_call = input$geno_call_col
      )
    })

    # Display KASP plate status
    output$kasp_color_code_stat <- DT::renderDT({
      req(color_coded_stat())
      DT::datatable(color_coded_stat(),
                    options = list(
                      pageLength = 10,
                      scrollX = TRUE
                    )
      )
    })

    return(color_coded) # return value to be used by next parameter.
  })
}

## To be copied in the UI
# mod_mv_kasp_color_ui("mv_kasp_color_1")

## To be copied in the server
# mod_mv_kasp_color_server("mv_kasp_color_1")
