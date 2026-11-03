#' Predictive Summary Statistics for KASP Genotyping Controls UI Function
#'
#' @description A Shiny module for analyzing and visualizing predictive statistics
#'   for positive controls in KASP genotyping data
#'
#' @param id Internal ID parameter for {shiny}
#'
#' @return A tagList containing the UI elements
#'
#' @noRd
#'
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel selectInput textInput actionButton icon div conditionalPanel column fluidRow numericInput plotOutput downloadButton
#' @importFrom bslib input_switch accordion accordion_panel card card_footer
#'
mod_mv_pred_sum_stat_ui <- function(id) {
  ns <- NS(id)
  ui_sidebar_layout(

    # SIDEBAR ---------------------------------------------------------------

    sidebar = tagList(
      selectInput(
        inputId = ns("snp_id"),
        label = "SNP ID Column",
        choices = NULL,
        multiple = FALSE
      ),
      selectInput(
        inputId = ns("geno_call_id"),
        label = "Genotype Call Column",
        choices = NULL,
        multiple = FALSE
      ),
      selectInput(
        inputId = ns("group_id"),
        label = label_help(
          "Positive Control Column",
          "Column with the expected genotype of each positive control. Predictions are scored against it."
        ),
        choices = NULL,
        multiple = FALSE
      ),
      bslib::input_switch(
        id = ns("rate_out_id"),
        label = label_help(
          "Show as proportions",
          "On: share of controls in each class. Off: raw counts."
        ),
        value = TRUE
      ),

      # Labels used in the call column
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          title = "Call labels",
          icon = icon("tags"),
          textInput(
            inputId = ns("group_unknown"),
            label = label_help(
              "Unverified Status Label",
              "Label used in the positive control column when the expected genotype is not known."
            ),
            value = "?"
          ),
          textInput(
            inputId = ns("blank"),
            label = label_help("NTC / Blank Label", "Label of the no-template-control (blank) wells."),
            value = "NTC"
          )
        )
      )
    ),

    # MAIN PANEL ------------------------------------------------------------

    # Summary table card
    ui_result_card(
      "Prediction summary for positive controls",
      full_screen = TRUE,
      DT::DTOutput(outputId = ns("predict_sum"))
    ),

    # Summary plot card, with SNP picker, Appearance and download
    ui_result_card(
      "Prediction summary plot",
      full_screen = TRUE,
      actions = tagList(
        div(
          style = "min-width: 200px; margin-bottom: -1rem;",
          selectInput(inputId = ns("SNP_id"), label = NULL, choices = NULL, multiple = FALSE, width = "100%")
        ),
        plot_options_button(
          bslib::layout_columns(
            numericInput(ns("textsize_id"), "Font Size", value = 12, min = 1, max = 30, step = 1),
            numericInput(ns("alpha_id"), "Bar Opacity", value = 1, min = 0, max = 1, step = 0.05)
          ),
          tags$p(class = "fw-semibold small mb-1", "Bar colours"),
          bslib::layout_columns(
            colourpicker::colourInput(ns("col_false"), "False", value = "firebrick3", returnName = TRUE),
            colourpicker::colourInput(ns("col_true"), "True", value = "cornflowerblue", returnName = TRUE),
            colourpicker::colourInput(ns("col_unverified"), "Unverified", value = "beige", returnName = TRUE)
          )
        ),
        plot_download_button(
          ns, "download_pred_plot",
          file_name = "Predictive_plot",
          width = 8, height = 5
        )
      ),
      plotOutput(outputId = ns("pred_plot"), height = "450px")
    ),

    # Per-sample status card, with plate picker
    ui_result_card(
      "Prediction status per sample",
      full_screen = TRUE,
      actions = div(
        style = "min-width: 260px; margin-bottom: -1rem;",
        selectInput(inputId = ns("pred_stat_dd"), label = NULL, choices = NULL, multiple = FALSE, width = "100%")
      ),
      DT::DTOutput(outputId = ns("predict_stat"))
    )
  )
}


#' Predictive Summary Statistics Server Functions
#'
#' @description Server logic for analyzing and visualizing predictive statistics
#'   for positive controls in KASP genotyping data
#'
#' @param id An ID string that corresponds with the ID used to call the module's UI function
#' @param kasp_data A reactive expression that returns the KASP data frame
#' @param color_code_res A reactive expression that returns the color-coded data
#'
#' @return A reactive expression with selected inputs for further processing
#'
#' @noRd
#'
#' @importFrom shiny moduleServer reactive observeEvent req validate need updateSelectInput renderPlot downloadHandler
#'
mod_mv_pred_sum_stat_server <- function(id, color_code_res, kasp_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # COLUMN DROPDOWNS ------------------------------------------------------

    # Column names of the KASP data
    drop_down <- reactive({
      req(kasp_data())
      colnames(kasp_data())
    })

    # Fill the dropdowns and guess the right column for each
    observeEvent(drop_down(), {
      req(drop_down())
      updateSelectInput(session,
                        inputId = "snp_id",
                        choices = drop_down(),
                        selected = safe_grep_match(pattern = 'SNP',
                                                   choices = drop_down())
      )

      updateSelectInput(session,
                        inputId = "geno_call_id",
                        choices = drop_down(),
                        selected = safe_grep_match(pattern = 'Call',
                                                   choices = drop_down()))

      updateSelectInput(session,
                        inputId = "group_id",
                        choices = drop_down(),
                        selected = safe_grep_match(pattern = 'Group',
                                                   choices = drop_down())
      )
    })

    # PREDICTION SUMMARY ----------------------------------------------------

    pred_sum_result <- reactiveVal()

    # Re-run pred_summary() whenever the inputs change
    observe({
      req(
        input$geno_call_id, input$group_id,
        input$snp_id, input$blank,
        input$group_unknown, color_code_res()
      )

      cc_res <- color_code_res()
      req(length(cc_res) >= 1)

      # Only go on if the chosen columns are in the data
      current_cols <- colnames(cc_res[[1]])
      req(
        input$snp_id %in% current_cols,
        input$group_id %in% current_cols,
        input$geno_call_id %in% current_cols
      )

      pred_result <- pred_summary(
        x = color_code_res(),
        Group_unknown = input$group_unknown,
        blank = input$blank,
        snp_id = input$snp_id,
        Group_id = input$group_id,
        geno_call = input$geno_call_id,
        rate_out = if(input$rate_out_id) TRUE else FALSE
      )

      pred_sum_result(pred_result)
    })

    # Summary table
    output$predict_sum <- DT::renderDT({
      validate(need(kasp_data(), msg_no_data))
      validate(need(if_ready(color_code_res()), msg_no_color))
      req(pred_sum_result())
      tryCatch(
        {
          req(pred_sum_result()$summ)
          DT::datatable(pred_sum_result()$summ, options = list(pageLength = 10, scrollX = TRUE))
        },
        error = function(e) {
          shinyWidgets::show_toast(
            title = "",
            text = paste("Couldn't show the prediction summary:", e$message),
            type = "error",
            timer = 2000,
            position = "bottom-end"
          )
          return(NULL)
        }
      )
    })

    # PER-SAMPLE STATUS -----------------------------------------------------

    # Plate names
    plates_field <- reactive({
      req(pred_sum_result())
      names(pred_sum_result()$plates)
    })

    # Fill the plate picker
    observe({
      req(plates_field())
      updateSelectInput(session, inputId = "pred_stat_dd", choices = plates_field())
    })

    # Status table for the selected plate
    output$predict_stat <- DT::renderDT({
      validate(need(kasp_data(), msg_no_data))
      validate(need(if_ready(color_code_res()), msg_no_color))
      req(pred_sum_result(), input$pred_stat_dd)
      tryCatch(
        DT::datatable(pred_sum_result()$plates[[input$pred_stat_dd]],
                      options = list(pageLength = 10, scrollX = TRUE)),
        error = function(e) NULL
      )
    })

    # SUMMARY PLOT ----------------------------------------------------------

    # One bar plot per SNP
    pred_sum_plot <- reactive({
      req(
        pred_sum_result(),
        input$col_false, input$col_true, input$col_unverified,
        input$textsize_id,
        input$alpha_id
      )

      tryCatch(
        {
          pred_summary_plot(
            x = as.data.frame(pred_sum_result()$summ),
            pred_cols = c(false = input$col_false,
                          true = input$col_true,
                          unverified = input$col_unverified),
            alpha = input$alpha_id,
            text_size = input$textsize_id
          )
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Couldn't draw the prediction plots",
            text = e$message,
            type = "error",
            showCloseButton = TRUE
          )
          return(NULL)
        }
      )
    })

    # SNP names from the plots
    name_snpid <- reactive({
      req(pred_sum_plot())
      tryCatch(
        {
          names(pred_sum_plot())
        },
        error = function(e) {
          shinyWidgets::show_toast(
            title = "",
            text = paste("Couldn't read the SNP names:", e$message),
            type = "error",
            timer = 2000,
            position = "bottom-end"
          )
        }
      )
    })

    # Fill the SNP picker
    observe({
      req(name_snpid())
      updateSelectInput(session, inputId = "SNP_id", choices = name_snpid())
    })

    # Plot for the selected SNP
    output$pred_plot <- renderPlot({
      validate(need(kasp_data(), msg_no_data))
      validate(need(if_ready(color_code_res()), msg_no_color))
      req(input$SNP_id, pred_sum_plot())
      tryCatch(
        print(pred_sum_plot()[[input$SNP_id]]),
        error = function(e) {
          shinyWidgets::show_toast(
            title = "",
            text = "Pick three colours (False, True and Unverified) under Appearance.",
            type = "error",
            timer = 2000
          )
        }
      )
    })

    # Toast once each time the plots are rebuilt
    observeEvent(pred_sum_plot(), {
      shinyWidgets::show_toast(
        title = "",
        text = "Plot rendered successfully",
        type = "success",
        timer = 2000
      )
    })

    # DOWNLOADS -------------------------------------------------------------

    # Plot for the selected SNP
    output$download_pred_plot <- downloadHandler(
      filename = function() {
        paste0(input$download_pred_plot_name, ".", input$download_pred_plot_type)
      },
      content = function(file) {
        req(pred_sum_plot(), input$SNP_id)

        plot_obj <- pred_sum_plot()[[input$SNP_id]]

        if (is.null(plot_obj)) {
          shinyWidgets::show_alert("Plot not available", "No prediction plot for this SNP. Pick another SNP from the list.", type = "warning")
          return()
        }

        ggplot2::ggsave(
          file,
          plot = plot_obj,
          device = input$download_pred_plot_type,
          width = as.numeric(input$download_pred_plot_width),
          height = as.numeric(input$download_pred_plot_height)
        )
      }
    )
  })
}

## To be copied in the UI
# mod_mv_pred_sum_stat_ui("mv_pred_sum_stat_1")

## To be copied in the server
# mod_mv_pred_sum_stat_server("mv_pred_sum_stat_1")
