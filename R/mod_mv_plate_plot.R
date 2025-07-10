#' Genotyping Plate Layout Visualization UI Function
#'
#' @description A Shiny module for visualizing KASP genotyping plate layouts
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel selectInput actionButton plotOutput div downloadButton textInput
#' @importFrom bslib accordion accordion_panel card card_footer
mod_mv_plate_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = ns("well_id"),
          label = "Select Column for Genotyping Plate wells",
          choices = NULL,
          multiple = FALSE,
          width = "100%"
        ),
        selectInput(
          inputId = ns("Geno_call_id"),
          label = "Select Column for Genotypes",
          choices = NULL,
          multiple = FALSE,
          width = "100%"
        ),
        selectInput(
          inputId = ns("snp_idd"),
          label = "Select Column for SNP ID",
          choices = NULL,
          multiple = FALSE,
          width = "100%"
        ),
        selectInput(
          inputId = ns("cplate_well"),
          choices = NULL,
          label = "Select Plate",
          multiple = FALSE,
          width = "100%"
        )
        # div(
        #   style = "display: flex; justify-content: center;",
        #   actionButton(
        #     inputId = ns("plot_btn"), label = "Generate Layout",
        #     class = "btn-primary", width = "70%"
        #   )
        # )
      ),
      mainPanel(
        bslib::accordion(
          open = "KASP Genotyping Plate Layout",
          bslib::accordion_panel(
            title = "KASP Genotyping Plate Layout",
            plotOutput(outputId = ns("plate_layout_plot"),width = '100%' ,height = '500px' ),
            bslib::card(card_footer(
              fluidRow(
                column(
                  3,
                  textInput(
                    inputId = ns("file_name"),
                    label = "Enter Filename",
                    value = "Plate layout 1"
                  )
                ),
                column(
                  3,
                  numericInput(
                    inputId = ns("width"),
                    label = "Set Plot Width",
                    value = 8, min = 1
                  )
                ), column(
                  3,
                  numericInput(
                    inputId = ns("height"),
                    label = "Set Plot Height",
                    value = 5, min = 1
                  )
                )
              ),
              downloadButton(
                outputId = ns("download_plateplot"),
                label = "Download Plot", class = "btn-success"
              )
            ))
          )
        )
      )
    )
  )
}

#' Genotyping Plate Layout Visualization Server Function
#'
#' @description Server logic for visualizing KASP genotyping plate layouts
#'
#' @param id An ID string that corresponds with the ID used to call the module's UI function
#' @param kasp_data A reactive expression that returns the KASP data frame
#' @param color_coded A reactive expression that returns the color-coded plate data
#'
#' @return None
#'
#' @noRd
#'
#' @importFrom shiny moduleServer reactive observeEvent req updateSelectInput renderPlot downloadHandler
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#'
mod_mv_plate_plot_server <- function(id, kasp_data, color_coded) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    sub_names <- reactive({
      req(kasp_data)
      colnames(kasp_data)
    })

    # Populate choices in subset
    observeEvent(req(sub_names()), {
      updateSelectInput(session, inputId = "snp_idd", choices = sub_names(), selected = "SNPID")
      # Plate layout
      updateSelectInput(session, inputId = "Geno_call_id", choices = sub_names(), selected = "Call")
      updateSelectInput(session, inputId = "well_id", choices = sub_names(), selected = "MasterWell")
    })

    names_col <- reactive({
      req(color_coded())
      names(color_coded())
    })

    # Populate dropdown for subplates
    observeEvent(names_col(), {
      req(names_col())
      updateSelectInput(session, inputId = "cplate_well", choices = names_col())
    })


    # Plate layout plot script.
    # Reactive for generating plate plot with error handling
    plot_plate_result <- reactive({
      tryCatch(
        {
          # Validate inputs
          req(color_coded(), input$snp_idd, input$Geno_call_id, input$well_id, input$cplate_well)

          # Generate plot
          plot_plate(
            x = color_coded()[input$cplate_well],
            well = input$well_id,
            geno_call = input$Geno_call_id,
            snp_id = input$snp_idd,
            filename = NULL,
            pdf = FALSE
          )
        },
        error = function(e) {
          show_alert(
            title = "Error!",
            text = paste("Failed to generate plot:", e$message),
            type = "error",
            showCloseButton = TRUE,
            timer = 5000
          )
          return(NULL)
        }
      )
    })

    # Display the plate layout plot
    output$plate_layout_plot <- renderPlot({
      req(plot_plate_result(), input$cplate_well)

      # Find matching plot (case-insensitive)
      plot_name <- grep(
        pattern = input$cplate_well,
        x = names(plot_plate_result()),
        ignore.case = TRUE,
        value = TRUE
      )

      if (length(plot_name) == 0) {
        showNotification("No plot found for selected plate", type = "error")
        return(NULL)
      }

      plot_obj <- plot_plate_result()[[plot_name]]

      if (is.null(plot_obj)) {
        showNotification("Plot object is NULL", type = "error")
        return(NULL)
      }

      # Show success message only after successful render
      shinyWidgets::show_toast(
        title = "Success",
        text = "Plot rendered successfully",
        type = "success",
        timer = 2000
      )


      print(plot_obj)
    })

    output$download_plateplot <- downloadHandler(
      filename = function() {
        req(input$file_name)
        paste0(input$file_name, ".pdf")  # Single PDF containing ALL plots
      },
      content = function(file) {
        req(plot_plate_result())  # Ensure plots exist

        tryCatch(
          {
            # Start PDF (onefile=TRUE ensures multi-page)
            pdf(file,
                width = input$width,
                height = input$height,
                onefile = TRUE)

            # Print ALL plots regardless of type
            for (plot_obj in plot_plate_result()) {
              if (inherits(plot_obj, "grob")) {
                grid::grid.draw(plot_obj)  # Handle grid graphics
              } else {
                print(plot_obj)  # Handle ggplot2/base R plots
              }
            }

            dev.off()
          },
          error = function(e) {
            showNotification(
              paste("Download failed:", e$message),
              type = "error",
              duration = 5
            )
          }
        )
      }
    )

    # # Download plate layout plot
    # output$download_plateplot <- downloadHandler(
    #   filename = function() {
    #     req(input$file_name)
    #     paste0(input$file_name, ".pdf")
    #   },
    #   content = function(file) {
    #     req(plot_plate_result())
    #
    #     tryCatch(
    #       {
    #         # Determine if we need grid.draw or print
    #         plot_obj <- plot_plate_result()[[1]]
    #
    #         if (inherits(plot_obj, "grob")) {
    #           pdf(file, width = input$width, height = input$height)
    #           grid::grid.draw(plot_obj)
    #           dev.off()
    #         } else if (inherits(plot_obj, "ggplot")) {
    #           ggsave(
    #             file,
    #             plot = plot_obj,
    #             width = input$width,
    #             height = input$height,
    #             units = "in"
    #           )
    #         } else {
    #           pdf(file, width = input$width, height = input$height)
    #           print(plot_obj)
    #           dev.off()
    #         }
    #       },
    #       error = function(e) {
    #         show_alert(
    #           title = "Download Error!",
    #           text = paste("Failed to save plot:", e$message),
    #           type = "error",
    #           showCloseButton = TRUE,
    #           timer = 5000
    #         )
    #       }
    #     )
    #   }
    # )
  })
}

## To be copied in the UI
# mod_mv_plate_plot_ui("mv_plate_plot_1")

## To be copied in the server
# mod_mv_plate_plot_server("mv_plate_plot_1")
