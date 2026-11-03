#' Marker Call Rate UI Function
#'
#' @description A Shiny module that summarizes and plots KASP marker call rates
#'   using \code{kasp_call_rate()}.
#'
#' @param id Internal ID parameter for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel selectizeInput textInput numericInput plotOutput downloadButton fluidRow column selectInput helpText
#' @importFrom bslib navset_card_tab nav_panel card_footer
#'
mod_mv_call_rate_ui <- function(id) {
  ns <- NS(id)
  ui_sidebar_layout(

    # SIDEBAR ---------------------------------------------------------------

    sidebar = tagList(
      selectizeInput(
        inputId = ns("snp_id"),
        label = "SNP ID Column",
        choices = NULL,
        multiple = FALSE
      ),
      selectizeInput(
        inputId = ns("plate_id"),
        label = "Master Plate ID Column",
        choices = NULL,
        multiple = FALSE
      ),
      selectizeInput(
        inputId = ns("geno_call_id"),
        label = "Genotype Call Column",
        choices = NULL,
        multiple = FALSE
      ),
      numericInput(
        inputId = ns("threshold"),
        label = label_help(
          "Call Rate Threshold (%)",
          "Markers whose call rate is below this value are flagged as failing."
        ),
        value = 95, min = 0, max = 100, step = 1
      ),
      # Labels used in the call column
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          title = "Call labels",
          icon = icon("tags"),
          textInput(
            inputId = ns("sep"),
            label = label_help(
              "Allele Separator",
              "Character between the two alleles in a call, e.g. ':' in A:G."
            ),
            value = ":"
          ),
          textInput(
            inputId = ns("no_call"),
            label = label_help("No Call Label", "Label used for wells that gave no call."),
            value = "?"
          ),
          textInput(
            inputId = ns("uncallable"),
            label = label_help("Uncallable Label", "Label used for wells the software could not call."),
            value = "Uncallable"
          ),
          textInput(
            inputId = ns("blank"),
            label = label_help(
              "NTC / Blank Label",
              "Label of the no-template-control (blank) wells. These are left out of call rates."
            ),
            value = "NTC"
          )
        )
      )
    ),

    # MAIN PANEL ------------------------------------------------------------

    # Headline numbers
    bslib::layout_column_wrap(
      width = "180px",
      fill = FALSE,
      bslib::value_box(
        title = "Markers",
        value = textOutput(ns("vb_markers"), inline = TRUE),
        showcase_layout = bslib::showcase_top_right(),
        showcase = icon("dna")
      ),
      bslib::value_box(
        title = "Plates",
        value = textOutput(ns("vb_plates"), inline = TRUE),
        showcase_layout = bslib::showcase_top_right(),
        showcase = icon("table-cells")
      ),
      bslib::value_box(
        title = "Overall call rate",
        value = textOutput(ns("vb_rate"), inline = TRUE),
        showcase_layout = bslib::showcase_top_right(),
        showcase = icon("percent")
      ),
      bslib::value_box(
        title = "Markers passing",
        value = textOutput(ns("vb_pass"), inline = TRUE),
        showcase_layout = bslib::showcase_top_right(),
        showcase = icon("circle-check")
      )
    ),
    # Plot and table tabs, with download on the right
    bslib::navset_card_tab(
      id = ns("call_rate_tabs"),
      full_screen = TRUE,
      bslib::nav_panel(
        title = "Marker Call Rate",
        value = "marker_call_rate",
        plotOutput(ns("marker_plot"), height = "450px")
      ),
      bslib::nav_panel(
        title = "No-call Reasons",
        value = "no_call_reason",
        plotOutput(ns("no_call_plot"), height = "450px")
      ),
      bslib::nav_panel(
        title = "Per Plate",
        value = "plate_call_rate",
        plotOutput(ns("plate_plot"), height = "450px")
      ),
      bslib::nav_panel(
        title = "Summary Table",
        value = "table",
        selectInput(
          inputId = ns("table_level"),
          label = NULL,
          choices = c("By marker" = "marker_summary",
                      "By marker and plate" = "plate_summary")
        ),
        DT::DTOutput(ns("call_rate_table"))
      ),
      bslib::nav_spacer(),
      bslib::nav_item(
        plot_download_button(
          ns, "download_plot",
          file_name = "call_rate", width = 8, height = 5,
          note = "Saves the plot in the open tab (the marker plot when the table is open)."
        )
      )
    )
  )
}


#' Marker Call Rate Server Function
#'
#' @param id An ID string that corresponds with the ID used to call the module's UI function.
#' @param kasp_data A reactive returning a data frame of KASP genotype calls.
#'
#' @noRd
#'
#' @importFrom shiny moduleServer reactive observe req updateSelectizeInput renderPlot downloadHandler
#'
mod_mv_call_rate_server <- function(id, kasp_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # COLUMN DROPDOWNS ------------------------------------------------------

    # Column names of the KASP data
    col_names <- reactive({
      req(kasp_data())
      colnames(kasp_data())
    })

    # Fill the dropdowns and guess the right column for each
    observe({
      req(col_names())
      updateSelectizeInput(session, inputId = "snp_id", server = TRUE,
                           choices = col_names(),
                           selected = safe_grep_match(pattern = "snp", choices = col_names()))
      updateSelectizeInput(session, inputId = "plate_id", server = TRUE,
                           choices = col_names(),
                           selected = safe_grep_match(pattern = "masterplate", choices = col_names()))
      updateSelectizeInput(session, inputId = "geno_call_id", server = TRUE,
                           choices = col_names(),
                           selected = safe_grep_match(pattern = "call", choices = col_names()))
    })

    # CALL RATES ------------------------------------------------------------

    # Run kasp_call_rate() whenever the inputs change
    call_rate_res <- reactive({
      req(kasp_data(), input$snp_id, input$plate_id, input$geno_call_id,
          input$sep, input$threshold)
      req(all(c(input$snp_id, input$plate_id, input$geno_call_id) %in% col_names()))

      tryCatch(
        kasp_call_rate(
          x = kasp_data(),
          snp_id = input$snp_id,
          plate_id = input$plate_id,
          geno_call = input$geno_call_id,
          sep = input$sep,
          blank = input$blank,
          no_call = input$no_call,
          uncallable = input$uncallable,
          threshold = input$threshold / 100
        ),
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Couldn't compute call rates",
            text = e$message,
            type = "error",
            showCloseButton = TRUE
          )
          NULL
        }
      )
    })

    # PLOTS AND TABLE -------------------------------------------------------

    # Call rate per marker
    output$marker_plot <- renderPlot({
      validate(need(kasp_data(), msg_no_data))
      req(call_rate_res())
      call_rate_res()$marker_call_rate
    }, res = 110)

    # Why wells had no call
    output$no_call_plot <- renderPlot({
      validate(need(kasp_data(), msg_no_data))
      req(call_rate_res())
      call_rate_res()$no_call_reason
    }, res = 110)

    # Call rate per marker and plate
    output$plate_plot <- renderPlot({
      validate(need(kasp_data(), msg_no_data))
      req(call_rate_res())
      call_rate_res()$plate_call_rate
    }, res = 110)

    # Summary table, by marker or by marker and plate
    output$call_rate_table <- DT::renderDT({
      validate(need(kasp_data(), msg_no_data))
      req(call_rate_res(), input$table_level)
      df <- call_rate_res()[[input$table_level]]
      req(df)
      num_cols <- intersect(c("call_rate", "ci_lower", "ci_upper"), colnames(df))
      df[num_cols] <- lapply(df[num_cols], function(v) round(100 * v, 1))
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    # VALUE BOXES -----------------------------------------------------------

    # Shows "-" until there are results
    tile <- function(f) {
      renderText({
        res <- if_ready(call_rate_res())
        if (is.null(res)) "-" else f(res)
      })
    }
    output$vb_markers <- tile(function(res) nrow(res$marker_summary))
    output$vb_plates <- tile(function(res) length(unique(res$plate_summary[[2]])))
    output$vb_rate <- tile(function(res) {
      m <- res$marker_summary
      paste0(round(100 * sum(m$called) / sum(m$n_wells), 1), "%")
    })
    output$vb_pass <- tile(function(res) {
      m <- res$marker_summary
      paste(sum(m$pass), "of", nrow(m))
    })

    # DOWNLOADS -------------------------------------------------------------

    # Plot in the open tab (marker plot when the table is open)
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(input$download_plot_name, ".", input$download_plot_type)
      },
      content = function(file) {
        req(call_rate_res())
        tab <- input$call_rate_tabs
        if (is.null(tab) || tab == "table") tab <- "marker_call_rate"
        ggplot2::ggsave(file,
                        plot = call_rate_res()[[tab]],
                        width = input$download_plot_width,
                        height = input$download_plot_height,
                        device = input$download_plot_type)
      }
    )
  })
}

## To be copied in the UI
# mod_mv_call_rate_ui("mv_call_rate_1")

## To be copied in the server
# mod_mv_call_rate_server("mv_call_rate_1", kasp_data = kasp_data)
