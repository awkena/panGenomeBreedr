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
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel selectInput textInput
#'   actionButton icon div conditionalPanel column fluidRow numericInput plotOutput downloadButton
#' @importFrom bslib input_switch accordion accordion_panel card card_footer
#'
mod_mv_pred_sum_stat_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
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
          label = "Positive Control Column",
          choices = NULL,
          multiple = FALSE
        ),
        textInput(
          inputId = ns("group_unknown"),
          label = "Unverified Status Label",
          value = "?"
        ),
        textInput(
          inputId = ns("blank"),
          label = "NTC / Blank Label",
          value = "NTC"
        ),
        bslib::input_switch(
          id = ns("rate_out_id"),
          label = "Return Predictions As Proportion",
          value = TRUE
        )

      ),
      mainPanel(
        bslib::accordion(
          open = TRUE,
          bslib::accordion_panel(
            title = "Predictive Summary for Positive Controls in KASP Genotype Data",
            DT::DTOutput(outputId = ns("predict_sum"))
          ), style = "margin-bottom: 15px;"
        ),
        bslib::accordion(
          bslib::accordion_panel(
            title = "Predictive Summary Plot",
            fluidRow(
              column(
                width = 4,
                selectInput(inputId = ns("SNP_id"), label = "Target SNP ID", choices = NULL, multiple = FALSE)
              ),
              column(
                width = 4,
                numericInput(inputId = ns("alpha_id"), label = "Bar Opacity", value = 1, min = 0, max = 1, step = 0.05)
              ),
              column(
                width = 4,
                numericInput(inputId = ns("height_id"), label = "Plot Height", value = 5, min = 1, max = 30, step = 1)
              )
            ),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  inputId = ns("pred_col_id"),
                  label = "Positive Control Colors",
                  choices = grDevices::colors(),
                  selected = c("firebrick3", "cornflowerblue", "beige"),
                  multiple = TRUE
                ),
                helpText(
                  tags$em("Required Order: False, True, Unverified")
                )
              ),
              column(
                width = 4,
                numericInput(inputId = ns("textsize_id"), label = "Plot Font Size", value = 12, min = 1, max = 30, step = 1)
              ),
              column(
                width = 4,
                numericInput(inputId = ns("width_id"), label = "Plot Width", value = 8, min = 1, max = 30, step = 1)
              )
            ), hr(),
            plotOutput(outputId = ns("pred_plot"), height = "500px"),
            bslib::card(card_footer(
              textInput(inputId = ns("file_name"), label = "Enter Filename", value = "Predictive_plot"),
              downloadButton(outputId = ns("download_pred_plot"), label = "Download Plot", class = "btn-success")
            ))
          ),
          style = "margin-bottom: 15px;"
        ),
        bslib::accordion(
          bslib::accordion_panel(
            title = "Predictive Status for Positive Controls in KASP Genotype Data",
            selectInput(inputId = ns("pred_stat_dd"), label = "Select Plate ID", choices = NULL, multiple = FALSE),
            DT::DTOutput(outputId = ns("predict_stat"))
          ),
          style = "margin-bottom: 15px;"
        )
      )
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
    # Script to give the user drop down inputs for the select input tabs.
    drop_down <- reactive({
      req(kasp_data)
      colnames(kasp_data)
    })

    # Populate the select input values.
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

    pred_sum_result <- reactiveVal() # empty reactive value to hold result

    # Store result in the reactive value when "Run" button is clicked
    observe({
      req(
        input$geno_call_id, input$group_id,
        input$snp_id, input$blank,
        input$group_unknown, color_code_res()
      )

      cc_res <- color_code_res()
      req(length(cc_res) >= 1)

      # Get the actual column names from the current data
      current_cols <- colnames(cc_res[[1]])

      # Only proceed if selected inputs exist in the current data
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

    # Extract plate names (ensure pred_sum_result() has valid data)
    plates_field <- reactive({
      req(pred_sum_result())
      names(pred_sum_result()$plates)
    })

    # Update dropdown for plates selection
    observe({
      req(plates_field())

     updateSelectInput(session, inputId = "pred_stat_dd", choices = plates_field() )

    })

    # Improved handling of predictive status display
    observeEvent(input$pred_stat_dd,{
      req(pred_sum_result())

      # Display predictive status with error handling
      output$predict_stat <- DT::renderDT({
        tryCatch(
          {
            DT::datatable(pred_sum_result()$plates[[input$pred_stat_dd]], options = list(pageLength = 10, scrollX = TRUE))
          },
          error = function(e) {
           # print(paste("Error displaying prediction status:", e$message))
            return(NULL)
          }
        )
      })
    })

    # Display prediction summary with error handling
    output$predict_sum <- DT::renderDT({
      req(pred_sum_result())
      tryCatch(
        {
          req(pred_sum_result()$summ)
          DT::datatable(pred_sum_result()$summ, options = list(pageLength = 10, scrollX = TRUE))
        },
        error = function(e) {
          shinyWidgets::show_toast(
            title = "",
            text = paste("Error displaying prediction summary:", e$message),
            type = "error",
            timer = 2000,
            position = "bottom-end"
          )
          return(NULL)
        }
      )
    })

    pred_sum_plot <- reactive({
      # Validate inputs
      req(
        pred_sum_result(),
        input$pred_col_id,
        input$textsize_id,
        input$height_id,
        input$alpha_id,
        input$width_id,length(input$pred_col_id) == 3
      )

      tryCatch(
        {
          # Generate plot
           pred_summary_plot(
            x = as.data.frame(pred_sum_result()$summ),
            pred_cols = c( false = input$pred_col_id[1],
                           true = input$pred_col_id[2],
                          unverified = input$pred_col_id[3]
                          ),
            alpha = input$alpha_id ,
            width = input$width_id ,
            height = input$height_id ,
            text_size = input$textsize_id
          )
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Error!",
            text = paste("Error generating prediction plots:", e$message),
            type = "error",
            showCloseButton = TRUE
          )
          return(NULL)
        }
      )
    })

    # Get unique SNP names from plots
    name_snpid <- reactive({
      req(pred_sum_plot())
      tryCatch(
        {
          names(pred_sum_plot())
        },
        error = function(e) {
          shinyWidgets::show_toast(
            title = "Error",
            text = paste("Error extracting SNP names:", e$message),
            type = "error",
            timer = 2000,
            position = "bottom-end"
          )
        }
      )
    })

    # Update dropdown for SNP selection with error handling
    observe({
      req(name_snpid())

    updateSelectInput(session, inputId = "SNP_id", choices = name_snpid())

    })



    # Render predictive plot
    observe({
      req(input$SNP_id, pred_sum_plot())

      output$pred_plot <- renderPlot({
        tryCatch({
        print(pred_sum_plot()[[input$SNP_id]])
          # Show success message
      shinyWidgets::show_toast(
        title = "",
        text = "Plot rendered successfully",
        type = "success",
        timer = 2000
      )
        },error = function(e){
          shinyWidgets::show_toast(
            title = "",
            text =  'Select 3 Colors for Predictive Plot',
            type = "error",
            timer = 2000
          )
        })

      })

    })



    # Download Plot with error handling
    output$download_pred_plot <- downloadHandler(
      filename = function() {
        paste0(input$file_name, ".pdf")
      },
      content = function(file) {
        req(pred_sum_plot(), input$SNP_id)

        # Find exact match
        plot_obj <- pred_sum_plot()[[input$SNP_id]]

        if (is.null(plot_obj)) {
          shinyWidgets::show_alert("Error", "Plot not found for selected SNP", type = "error")
          return()
        }

      ggplot2::ggsave(
          file,
          plot = plot_obj,
          width = as.numeric(input$width_id),
          height = as.numeric(input$height_id)
        )
      }
    )


  })
}

## To be copied in the UI
# mod_mv_pred_sum_stat_ui("mv_pred_sum_stat_1")

## To be copied in the server
# mod_mv_pred_sum_stat_server("mv_pred_sum_stat_1")
