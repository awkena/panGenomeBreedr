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
#' @importFrom DT DTOutput
#' @importFrom grDevices colors
mod_mv_pred_sum_stat_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = ns("snp_id"),
          label = "Select Column for SNP ID",
          choices = NULL,
          multiple = FALSE
        ),
        selectInput(
          inputId = ns("geno_call_id"),
          label = "Select Column for Genotype Calls",
          choices = NULL,
          multiple = FALSE
        ),
        selectInput(
          inputId = ns("group_id"),
          label = "Select Column for Group ID",
          choices = NULL,
          multiple = FALSE
        ),
        textInput(
          inputId = ns("group_unknown"),
          label = "Unverified Genotype Indicator",
          value = "?"
        ),
        textInput(
          inputId = ns("blank"),
          label = "Specify No Template Control Calls",
          value = "NTC"
        ),
        bslib::input_switch(
          id = ns("rate_out_id"),
          label = "Return Proportions",
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
        shiny::conditionalPanel(
          condition = "input.rate_out_id == true",
          ns = ns,
          bslib::accordion(
            bslib::accordion_panel(
              title = "Predictive Summary Plot",
              fluidRow(
                column(
                  width = 4,
                  selectInput(inputId = ns("SNP_id"), label = "Choose SNP ID", choices = NULL, multiple = FALSE)
                ),
                column(
                  width = 4,
                  selectInput(
                    inputId = ns("pred_col_id"), label = "Pos. Ctrl Colors (F | T | U)",
                    choices = colors(),
                    selected = c("red", "blue", "orange"),
                    multiple = TRUE
                  )
                ),
                column(
                  width = 4,
                  numericInput(inputId = ns("alpha_id"), label = "Adjust Alpha Value", value = 1, min = 0, max = 1, step = 0.05)
                )
              ),
              fluidRow(
                column(
                  width = 4,
                  numericInput(inputId = ns("height_id"), label = "Adjust Plot Height", value = 5, min = 1, max = 30, step = 1)
                ),
                column(
                  width = 4,
                  numericInput(inputId = ns("width_id"), label = "Adjust Plot Width", value = 8, min = 1, max = 30, step = 1)
                ),
                column(
                  width = 4,
                  numericInput(inputId = ns("textsize_id"), label = "Text Size", value = 12, min = 1, max = 30, step = 1)
                )
              ),
              plotOutput(outputId = ns("pred_plot"), height = "500px"),
              bslib::card(card_footer(
                textInput(inputId = ns("file_name"), label = "Enter Filename", value = "Predictive_plot"),
                downloadButton(outputId = ns("download_pred_plot"), label = "Download Plot", class = "btn-success")
              ))
            ),
            style = "margin-bottom: 15px;"
          )
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
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom shinyjs delay
#' @importFrom shinyWidgets show_alert
#' @importFrom DT renderDT datatable
#' @importFrom grid grid.draw
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
      updateSelectInput(session, inputId = "snp_id", choices = drop_down(), selected = drop_down()[grep(pattern = "SNP", x = drop_down(), ignore.case = TRUE)])
      updateSelectInput(session, inputId = "geno_call_id", choices = drop_down(), selected = drop_down()[grep(pattern = "Call", x = drop_down(), ignore.case = TRUE)])
      updateSelectInput(session, inputId = "group_id", choices = drop_down(), selected = drop_down()[grep(pattern = "Group", x = drop_down(), ignore.case = TRUE)])
    })

    pred_sum_result <- reactiveVal()

    # Add processing state tracker
    #  processing <- reactiveVal(FALSE)

    # Store result in the reactive value when "Run" button is clicked
    observe({
      req(
        input$geno_call_id, input$group_id,
        input$snp_id, input$blank,
        input$group_unknown, color_code_res()
      )

      # Use shinyjs::delay to ensure UI updates before processing
      shinyjs::delay(100, {
        tryCatch(
          {
            # Validate required inputs
            validate(need(input$group_id, "Group ID is required"))
            validate(need(input$geno_call_id, "Genotype call ID is required"))
            validate(need(input$snp_id, "SNP ID is required"))
            validate(need(color_code_res(), "Color code results are required"))

            # Execute prediction function
            pred_result <- pred_summary(
              x = color_code_res(),
              Group_unknown = input$group_unknown,
              blank = input$blank,
              snp_id = input$snp_id,
              Group_id = input$group_id,
              geno_call = input$geno_call_id,
              rate_out = if (is.null(input$rate_out_id)) FALSE else input$rate_out_id
            )

            # Store result if valid
            if (!is.null(pred_result)) {
              pred_sum_result(pred_result)
            } else {
              show_alert(
                title = "Warning!",
                text = "Prediction analysis produced no results. Please check your inputs.",
                type = "warning",
                showCloseButton = TRUE,
                timer = 5000
              )
            }
          },
          error = function(e) {
            error_msg <- paste("Error in prediction analysis:", e$message)
            print(error_msg)
            show_alert(
              title = "Error!",
              text = error_msg,
              type = "error",
              showCloseButton = TRUE,
              timer = NULL
            )
          },
          warning = function(w) {
            warning_msg <- paste("Warning in prediction analysis:", w$message)
            print(warning_msg)
            show_alert(
              title = "Warning!",
              text = warning_msg,
              type = "warning",
              showCloseButton = TRUE,
              timer = 8000
            )
          }
        )
      })
    })

    # Extract plate names (ensure pred_sum_result() has valid data)
    plates_field <- reactive({
      req(pred_sum_result())
      names(pred_sum_result()$plates)
    })

    # Update dropdown for plates selection with error handling
    observeEvent(plates_field(), {
      tryCatch(
        {
          plates <- plates_field()
          if (length(plates) > 0) {
            updateSelectInput(session, inputId = "pred_stat_dd", choices = plates, selected = plates[1])
          } else {
            updateSelectInput(session, inputId = "pred_stat_dd", choices = character(0))
            show_alert(
              title = "Warning",
              text = "No plates found in the prediction results.",
              type = "warning",
              showCloseButton = TRUE
            )
          }
        },
        error = function(e) {
          print(paste("Error updating plates dropdown:", e$message))
        }
      )
    })

    # Get predictive status based on user's selected plate
    pred_stat_prime <- reactive({
      req(pred_sum_result(), input$pred_stat_dd)

      # Ensure selected plate exists before accessing
      if (input$pred_stat_dd %in% names(pred_sum_result()$plates)) {
        return(pred_sum_result()$plates[[input$pred_stat_dd]])
      } else {
        print("Error: Selected plate not found in pred_sum_result")
        return(NULL)
      }
    })

    # Improved handling of predictive status display
    observe({
      req(pred_stat_prime())

      # Display predictive status with error handling
      output$predict_stat <- DT::renderDT({
        req(pred_stat_prime())
        tryCatch(
          {
            DT::datatable(pred_stat_prime(), options = list(pageLength = 10, scrollX = TRUE))
          },
          error = function(e) {
            print(paste("Error displaying prediction status:", e$message))
            return(NULL)
          }
        )
      })

      # Show success message with appropriate delay
        shinyWidgets::show_toast(
          text = "Predictive Status Generated Successfully",
          title = "Success!",
          type = 'success',
          timer = 2000
        )

    })

    # Display prediction summary with error handling
    output$predict_sum <- DT::renderDT({
      req(pred_sum_result())
      tryCatch(
        {
          req(pred_sum_result()$summ)
          DT::datatable(as.data.frame(pred_sum_result()$summ), options = list(pageLength = 10, scrollX = TRUE))
        },
        error = function(e) {
          print(paste("Error displaying prediction summary:", e$message))
          return(NULL)
        }
      )
    })

    pred_sum_plot <- reactive({
      req(pred_sum_result())
      tryCatch(
        {
          # Validate inputs
          req(
            pred_sum_result()$summ,
            input$pred_col_id,
            input$textsize_id,
            input$height_id,
            input$alpha_id,
            input$width_id
          )

          validate(
            need(length(input$pred_col_id) == 3,
                 "Please select exactly 3 columns for predictions")
          )

          # Generate plot
           pred_summary_plot(
            x = as.data.frame(pred_sum_result()$summ),
            pred_cols = stats::setNames(input$pred_col_id[1:3], c('false', 'true', 'unverified')),
            alpha = input$alpha_id ,
            width = input$width_id ,
            height = input$height_id ,
            text_size = input$textsize_id
          )
        },
        error = function(e) {
          error_msg <- paste("Error generating prediction plots:", e$message)
          message(error_msg)

          show_alert(
            title = "Error!",
            text = error_msg,
            type = "error",
            showCloseButton = TRUE
          )
          return(NULL)
        }
      )
    })
    # Get unique SNP names from summary with error handling
    name_snpid <- reactive({
      req(pred_sum_plot())
      tryCatch(
        {
          names(pred_sum_plot())
        },
        error = function(e) {
          print(paste("Error extracting SNP names:", e$message))
          return(character(0))
        }
      )
    })

    # Update dropdown for SNP selection with error handling
    observe({
      req(pred_sum_plot())
      tryCatch(
        {
          snp_names <- name_snpid()
          if (length(snp_names) > 0) {
            updateSelectInput(session, inputId = "SNP_id", choices = snp_names)
          } else {
            updateSelectInput(session, inputId = "SNP_id", choices = character(0))
            show_alert(
              title = "Warning",
              text = "No SNP IDs found in the prediction results.",
              type = "warning",
              showCloseButton = TRUE
            )
          }
        },
        error = function(e) {
          print(paste("Error updating SNP dropdown:", e$message))
        }
      )
    })



    # # Dynamic plots with error handling
    # plot_me <- reactive({
    #   req(pred_sum_plot(), input$SNP_id)
    #   tryCatch(
    #     {
    #       if (input$SNP_id %in% names(pred_sum_plot())) {
    #         pred_sum_plot()[input$SNP_id]
    #       } else {
    #         print("Error: Selected SNP not found in pred_sum_plot")
    #         return(NULL)
    #       }
    #     },
    #     error = function(e) {
    #       print(paste("Error retrieving plot:", e$message))
    #       return(NULL)
    #     }
    #   )
    # })

    # Render predictive plot
    observe({
      req(input$SNP_id, pred_sum_plot())

      output$pred_plot <- renderPlot({
        # check and pull exaact name from list plot names
        plot_name <- grep(
          pattern = input$SNP_id,
          x = names(pred_sum_plot()),
          ignore.case = TRUE,
          value = TRUE
        )

        # Validate
        if (length(plot_name) == 0) {
          showNotification("No plot found for selected SNP", type = "error")
          return(NULL)
        }

        # Get the plot object
        plot_obj <- pred_sum_plot()[[plot_name]]

        # Validate plot object
        if (is.null(plot_obj)) {
          showNotification("Plot object is NULL", type = "error")
          return(NULL)
        }

        tryCatch(
          {
            # Render the plot
            print(plot_obj)

            # Show success message
            shinyWidgets::show_toast(
              title = "Success",
              text = "Plot rendered successfully",
              type = "success",
              timer = 2000
            )
          },
          error = function(e) {
            error_msg <- paste("Plot rendering failed:", e$message)
            cat(error_msg, "\n")
            showNotification(error_msg, type = "error")
            return(NULL)
          }
        )
      })
    }) |>
      bindEvent(input$SNP_id)  # only change when input for snpid changes


    # Download Plot with error handling
    output$download_pred_plot <- downloadHandler(
      filename = function() {
        paste0(input$file_name, ".pdf")
      },
      content = function(file) {
        req(pred_sum_plot(), input$SNP_id)

        # Find exact match
        plot_obj <- pred_sum_plot()[[grep(input$SNP_id, names(pred_sum_plot()), ignore.case = TRUE, value = TRUE)]]

        if (is.null(plot_obj)) {
          show_alert("Error", "Plot not found for selected SNP", type = "error")
          return()
        }

        ggsave(
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
