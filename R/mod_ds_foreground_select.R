#' ds_foreground_select UI Function
#'
#' @description A shiny Module for selecting foreground genotypes based on marker info.
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
        title = "Foreground select",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(width = '350px',
         bslib::card(
           bslib::card_header('Set Parameters'),
           bslib::card_body(
             fileInput(inputId = ns("geno_data"),
                       label = "Upload Genotypic Data (geno)",
                       accept = c('.xlsx' ,'.xls',".csv")),

             fileInput(inputId = ns("marker_info"),
                       label =  "Upload Marker Info File",
                       accept = c('.xlsx' ,'.xls',".csv")
             ),

             textInput(inputId = ns("sep"),
                       label = "Allele Separator",
                       value = ":"),

             selectInput(inputId = ns("fore_marker_col"),
                         label = "Foreground Marker Column",
                         choices = NULL),

             selectInput(inputId = ns("fav_allele_col"),
                         label = "Favorable Allele Column",
                         choices = NULL),

             selectInput(inputId = ns("alt_allele_col"),
                         label = "Alternative Allele Column",
                         choices = NULL),

             selectInput(inputId = ns("select_type"),
                         label = "Selection Type",
                         choices = c("homo", "hetero", "both"),
                         selected = "homo")
           ),
           bslib::card_footer(
            actionButton(inputId = ns("run_analysis"),
                         label = "Run Foreground Selection",
                         width = "100%",
                         icon = icon("rocket"),
                        # class = 'btn-info'
                        style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
                        `onmouseover` = "this.style.backgroundColor='#145214'",
                        `onmouseout` = "this.style.backgroundColor='forestgreen'"
                         )
            )
         )
          ),

          # Main panel.
          bslib::card(max_height = '700px',
            bslib::card_header("Configure UpSet Plot", class = "bg-primary text-white"),
            bslib::card_body(
              fluidRow(
                column(3,
                       textInput(
                         inputId = ns("mainbar_y_label"),
                         label = "Main Bar Y-axis Label",
                         value = "Locus Intersection Size"
                       )
                ),
                column(3,
                       textInput(
                         inputId = ns("sets_x_label"),
                         label = "Sets X-axis Label",
                         value = "Locus Size"
                       )
                ),
                column(3,
                       numericInput(
                         inputId = ns("text_scale"),
                         label = "Text Scale",
                         value = 1.2,
                         min = 0.5,
                         max = 3,
                         step = 0.1
                       )
                ),
                column(3,
                       selectInput(
                         inputId = ns("plot_type"),
                         label = "Metadata Plot Type",
                         choices = c("text"),
                         selected = "text"
                       )
                )
              ),
              fluidRow(
                column(3,
                       numericInput(
                         inputId = ns("plot_assign"),
                         label = "Assign (Position)",
                         value = 8,
                         min = 1,
                         step = 1
                       )
                ),
                column(3,
                       textInput(
                         inputId = ns("plot_column"),
                         label = "Metadata Column",
                         value = "locus"
                       )
                ),
                column(3,
                       selectInput(
                         inputId = ns("plot_colors"),
                         label = "Colors (comma-separated)",
                         selected = "darkblue",
                         choices = grDevices::colors()

                       )
                )
              ),

              bslib::input_switch(id = ns("configure"), label = "Findlines", value = FALSE),
              conditionalPanel(
                condition = paste0('input["', ns("configure"), '"] == true'),



                bslib::card(
                  class = "shadow p",
                  bslib::card_body(
                    fluidRow(

                      # Filter Controls Column
                      column(
                        width = 6,
                        tagList(
                          selectInput(
                            inputId = ns("present"),
                            label = "Filter SNPs to be present",
                            choices = NULL,
                            multiple = TRUE,
                            width = '100%'
                          ),
                          selectInput(
                            inputId = ns("absent"),
                            label = "Filter SNPs to be absent",
                            choices = NULL,
                            multiple = TRUE,
                            width = '100%'

                          ),
                          div(
                            style = "display: flex; justify-content: center;",
                            actionButton(
                              inputId = ns("Search"),
                              label = "Get lines",
                              icon = icon("search"),
                             # class = "btn",
                              width = "60%",
                             style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
                             `onmouseover` = "this.style.backgroundColor='#145214'",
                             `onmouseout` = "this.style.backgroundColor='forestgreen'"
                            )


                          )
                        )
                      ),

                      # Display Results Column
                      column(
                        width = 6,
                        bslib::card(
                          #bslib::card_header("Matching Genotypes"),
                          style = "height: 200px; overflow-y: auto;",
                          verbatimTextOutput(outputId = ns("display_result"))
                        )

                      )

                    )
                  )
                )

               )
            )
          ),
          # Results display section
          bslib::card(
            max_height = '800px',
            height = "auto",

            bslib::card_header(tags$b("Results")),
            bslib::card_body(
              navset_card_tab(
                id = ns("results_tabs"),

                # Data Table Tab
                bslib::nav_panel(
                  title = "Data Table",
                  icon = icon("table"),
                  DT::DTOutput(ns("result_table"), height = '800px')
                ),

                # Upset Plot Tab
                bslib::nav_panel(
                  title = "Upset Plot",
                  icon = icon("chart-line"),
                  # Upset Plot Output
                  plotOutput(ns("result_plot"), height = '800px')
                )
              )
            )
          )

        )
      )
    )
  )
}



#' ds_foreground_select Server Function
#'
#' @noRd
mod_ds_foreground_select_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Reactive values to hold uploaded data
    geno_data <- reactive({
      req(input$geno_data)

      read_mapfile(input$geno_data$datapath)
    })

    marker_info <- reactive({
      req(input$marker_info)

      read_mapfile(input$marker_info$datapath)
    })

    # Dynamically update select input choices when marker info is uploaded
    observeEvent(marker_info(), {
      cols <- colnames(marker_info())
      updateSelectInput(session, "fore_marker_col",
                        choices = cols,
                        selected = safe_grep_match(pattern = 'marker' ,choices = cols ))

      updateSelectInput(session, "fav_allele_col",
                        choices = cols,
                        selected = safe_grep_match(pattern = 'fav' ,choices = cols ))

      updateSelectInput(session,
                        "alt_allele_col",
                        choices = cols,
                        selected =  safe_grep_match(pattern = 'alt' ,choices = cols ) )
    })

    # Run analysis when button is clicked
    result_data <- eventReactive(input$run_analysis, {
      req(
        geno_data(), marker_info(),
        input$fore_marker_col, input$fav_allele_col,
        input$alt_allele_col, input$sep, input$select_type
      )

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Running foreground analysis... Please wait."
      )

      result <- tryCatch({
        foreground_select(
          geno_data = geno_data(),
          fore_marker_info = marker_info(),
          fore_marker_col = input$fore_marker_col,
          fav_allele_col = input$fav_allele_col,
          alt_allele_col = input$alt_allele_col,
          select_type = input$select_type,
          sep = check_sep(input$sep)
        )

        # shinyWidgets::show_toast(
        #   title = "Successful",
        #   text = '',
        #   type = "success"
        # )

      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Error",
          text = paste("An error occurred during foreground analysis:", e$message),
          type = "error"
        )
        return(NULL)
      }, finally = {
        shinyjs::delay(ms = 2000, {
          shinybusy::remove_modal_spinner()
        })
      })

      result
    })


    # Render DT Table
    output$result_table <- DT::renderDT({
      req(result_data())
      DT::datatable(result_data(), options = list(pageLength = 15))
    })


    # Generating the upset plot.
    upset_result <- reactive({
      req(
        result_data(), marker_info(),
        input$mainbar_y_label, input$sets_x_label,
        input$text_scale, input$plot_type,
        input$plot_assign, input$plot_column,
        input$plot_colors
      )

      tryCatch({
        run_upset_plot(
          foreground_matrix = result_data(),
          marker_info = marker_info(),
          mainbar_y_label = input$mainbar_y_label,
          sets_x_label = input$sets_x_label,
          text_scale = input$text_scale,
          plot_type = input$plot_type,
          assign = input$plot_assign,
          column = input$plot_column,
          colors = input$plot_colors
        )
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Error",
          text = paste("Upset plot failed:", e$message),
          type = "error"
        )
        return(NULL)
      })
    })

    # Render upset plot
    output$result_plot <- renderPlot({
      req(upset_result())
      print(upset_result())
    })

    # Show toast once plot is ready
    observeEvent(upset_result(), {
      shinyWidgets::show_toast(
        title = "Plot Ready",
        text = "Upset plot has been successfully rendered!",
        type = "success",
        position = "top-end",
        timer = 3000
      )
    })


    # Update select input
    observe({
      req(result_data())

      updateSelectInput(session ,inputId = 'present',
                        choices = colnames(result_data()),
                        selected = character(0))

      updateSelectInput(session ,inputId = 'absent',
                        choices = colnames(result_data()),
                        selected = character(0))
    })



    found_lines <- eventReactive(input$Search,{
      #input$present , input$absent,
      req(result_data())

      find_lines(mat = result_data(),present = input$present ,absent = input$absent)

    })



    # render result.
    output$display_result <- renderPrint({
      req(found_lines())
      paste0('line',found_lines())
    })

  })
}



## To be copied in the UI
# mod_ds_foreground_select_ui("ds_foreground_select_1")

## To be copied in the server
# mod_ds_foreground_select_server("ds_foreground_select_1")

