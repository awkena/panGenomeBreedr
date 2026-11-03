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
  ui_sidebar_layout(

    # SIDEBAR ---------------------------------------------------------------

    sidebar = tagList(
      selectInput(
        inputId = ns("Subset_names"),
        label = label_help(
          "Plate Subset Column",
          "Column that identifies each plate-by-marker subset (the 'plates' column added on upload). Calls are colour-coded per subset."
        ),
        choices = NULL, multiple = FALSE
      ),
      selectInput(
        inputId = ns("geno_call_col"),
        label = "Genotype Call Column",
        choices = NULL, multiple = FALSE
      ),
      textInput(
        inputId = ns("sep"),
        label = label_help(
          "Genotype Separator",
          "Character between the two alleles in a call, e.g. ':' in A:G."
        ),
        value = ":"
      ),
      bslib::accordion(
        open = FALSE,
        multiple = TRUE,

        # Labels used in the call column
        bslib::accordion_panel(
          title = "Call labels",
          icon = icon("tags"),
          textInput(
            inputId = ns("uncallable"),
            label = label_help("Uncallable Label", "Label used for wells the software could not call."),
            value = "Uncallable"
          ),
          textInput(
            inputId = ns("unused"),
            label = label_help("Unused Well Label", "Label used for empty or unused wells."),
            value = "?"
          ),
          textInput(
            inputId = ns("blank"),
            label = label_help("NTC / Blank Label", "Label of the no-template-control (blank) wells."),
            value = "NTC"
          ),
          textInput(
            inputId = ns("others"),
            label = label_help(
              "Non-genotype Labels",
              "Other labels in the call column that are not genotypes, separated by commas."
            ),
            value = "Missing , Bad , Dupe , Over , Short"
          )
        ),

        # FAM, HEX and het colours
        bslib::accordion_panel(
          title = "Genotype colours",
          icon = icon("palette"),
          colourpicker::colourInput(ns("col_fam"), "FAM homozygote", value = "blue", returnName = TRUE),
          colourpicker::colourInput(ns("col_hex"), "HEX homozygote", value = "gold", returnName = TRUE),
          colourpicker::colourInput(ns("col_het"), "Heterozygote", value = "forestgreen", returnName = TRUE)
        )
      )
    ),

    # MAIN PANEL ------------------------------------------------------------

    # Plate status card
    ui_result_card(
      "Plate status after colour coding",
      full_screen = TRUE,
      DT::DTOutput(outputId = ns("kasp_color_code_stat"))
    ),

    # Colour-coded calls card, with a plate picker in the header
    ui_result_card(
      "Colour-coded genotype calls",
      full_screen = TRUE,
      actions = div(
        style = "min-width: 260px; margin-bottom: -1rem;",
        selectInput(
          inputId = ns("subplates"),
          choices = NULL,
          label = NULL, multiple = FALSE, width = "100%"
        )
      ),
      DT::DTOutput(outputId = ns("kasp_color_code"))
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

    # COLUMN DROPDOWNS ------------------------------------------------------

    # Get colnames for drop down
    sub_names <- reactive({
      req(kasp_data())
      colnames(kasp_data())
    })

    # Fill the dropdowns and guess the right column for each
    observe({
      req(sub_names())
      updateSelectInput(session,
                        inputId = "Subset_names",
                        choices = sub_names(),
                        selected = safe_grep_match(pattern = 'plates',
                                                   choices = sub_names())
      )

      updateSelectInput(session,
                        inputId = "geno_call_col",
                        choices = sub_names(),
                        selected = safe_grep_match(pattern = 'call',
                                                   choices = sub_names())
      )
    })

    # COLOUR CODING ---------------------------------------------------------

    # Colour code the calls per plate
    color_coded <- reactive({
      req(
        kasp_data(),
        input$Subset_names,
        input$geno_call_col,
        input$col_fam, input$col_hex, input$col_het,
        input$uncallable,
        input$unused,
        input$blank,
        input$others,
        input$sep
      )

      tryCatch({
        kasp_color(
          x = kasp_data(),
          uncallable = trimws(input$uncallable),
          unused = trimws(input$unused),
          blank = trimws(input$blank),
          others = trimws(unlist(strsplit(input$others, split = ","))),
          sep = check_sep(input$sep),
          subset = input$Subset_names,
          geno_call = input$geno_call_col,
          assign_cols = c(
            FAM = input$col_fam,
            HEX = input$col_hex,
            het = input$col_het
          )
        )
      }, error = function(e) {
        shinyWidgets::show_alert("Couldn't colour-code the calls", paste0(e$message, " Check the genotype call column and separator."), type = "error")
        NULL
      })
    })

    # Toast each time the calls are re-coded
    observeEvent(color_coded(), {
      shinyWidgets::show_toast(
        title = "",
        type = "success",
        text = "Genotype calls colour-coded",
        timer = 2000,
        position = "bottom-end"
      )
    })

    # Plate names
    names_col <- reactive({
      req(color_coded())
      names(color_coded())
    })

    # Fill the plate picker
    observe({
      req(names_col())
      updateSelectInput(session, inputId = "subplates", choices = names_col())
    })

    # Colour-coded calls table
    output$kasp_color_code <- DT::renderDT({
      validate(need(kasp_data(), msg_no_data))
      validate(need(if_ready(color_coded()), "No colour-coded calls. Check the column and label settings."))
      req(input$subplates)
      DT::datatable(as.data.frame(color_coded()[[input$subplates]]),
                    options = list(pageLength = 10, scrollX = TRUE)
      )
    })

    # PLATE STATUS ----------------------------------------------------------

    # Status of plates
    color_coded_stat <- reactive({
      req(kasp_data(), input$Subset_names, input$geno_call_col)

      kasp_color_stat.give(
        x = kasp_data(),
        subset = input$Subset_names,
        geno_call = input$geno_call_col
      )
    })

    # Plate status table
    output$kasp_color_code_stat <- DT::renderDT({
      validate(need(kasp_data(), msg_no_data))
      req(color_coded_stat())
      DT::datatable(color_coded_stat(),
                    options = list(
                      pageLength = 10,
                      scrollX = TRUE
                    )
      )
    })

    # Colour-coded data goes to the QC plot and prediction modules
    return(color_coded)
  })
}

## To be copied in the UI
# mod_mv_kasp_color_ui("mv_kasp_color_1")

## To be copied in the server
# mod_mv_kasp_color_server("mv_kasp_color_1")
