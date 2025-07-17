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
                        selected = "homo"),

            actionButton(inputId = ns("run_analysis"),
                         label = "Run Foreground Selection",
                         icon = icon("rocket"),
                         class = 'btn-info')
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


              bslib::card(
           #     height = '500px',
                bslib::card_header(tags$b("Findlines")),
                bslib::card_body(
                  fluidRow(

                    # Filter Controls Column
                    column(
                      width = 5,
                      tagList(
                        selectInput(
                          inputId = ns("present"),
                          label = "Filter SNPs to be present",
                          choices = NULL,
                          multiple = TRUE
                        ),
                        selectInput(
                          inputId = ns("absent"),
                          label = "Filter SNPs to be absent",
                          choices = NULL,
                          multiple = TRUE
                        ),
                        actionButton(
                          inputId = ns("Search"),
                          icon = icon("search"),
                          label = "Search"
                        )
                      )
                    ),

                    # Display Results Column
                    column(
                      width = 7,
                      verbatimTextOutput(outputId = ns("display_result"))
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

      ext <- tools::file_ext(input$geno_data$name)
      if (ext == "csv") {
        read.csv(input$geno_data$datapath)|> as.data.frame()
      } else if (ext %in% c("xlsx", "xls")) {
        readxl::read_excel(input$geno_data$datapath)|> as.data.frame()
      }
    })

    marker_info <- reactive({
      req(input$marker_info)

      ext <- tools::file_ext(input$marker_info$name)

      if (ext == "csv") {
        read.csv(input$marker_info$datapath) |> as.data.frame()
      } else if (ext %in% c("xlsx", "xls")) {
        readxl::read_excel(input$marker_info$datapath,.name_repair = 'minimal') |> as.data.frame()
      }
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
      req(geno_data(), marker_info(),
          input$fore_marker_col, input$fav_allele_col, input$alt_allele_col,
          input$sep, input$select_type)

      foreground_select(
        geno_data = geno_data(),
        fore_marker_info = marker_info(),
        fore_marker_col = input$fore_marker_col,
        fav_allele_col = input$fav_allele_col,
        alt_allele_col = input$alt_allele_col,
        select_type = input$select_type,
        sep = check_sep(input$sep)
      )
    })

    # Render DT Table
    output$result_table <- DT::renderDT({
      req(result_data())
      DT::datatable(result_data(), options = list(pageLength = 15))
    })


    # Generating the upset plot.
    upset_result <- reactive({
      req(result_data() ,marker_info(),
          input$mainbar_y_label,input$sets_x_label,
          input$text_scale,input$plot_type, input$plot_assign,
          input$plot_column,input$plot_colors
          )

      run_upset_plot(foreground_matrix = result_data(),
                     marker_info = marker_info() ,
                     mainbar_y_label = input$mainbar_y_label,
                     sets_x_label = input$sets_x_label ,
                     text_scale = input$text_scale,
                     plot_type = input$plot_type ,
                     assign = input$plot_assign ,
                     column = input$plot_column ,
                     colors =  input$plot_colors)
    })

    # Render upset plot
    output$result_plot <- renderPlot({
      req(upset_result())
       print(upset_result())
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

    observe({
      req(found_lines())
      print(found_lines())
    })

    # # render result.
    # output$result_display <- renderText({
    #   req(found_lines())
    #   print(found_lines())
    # })

  })
}



## To be copied in the UI
# mod_ds_foreground_select_ui("ds_foreground_select_1")

## To be copied in the server
# mod_ds_foreground_select_server("ds_foreground_select_1")

