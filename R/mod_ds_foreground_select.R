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
          sidebar = bslib::sidebar(
            width = 380,
            position = "left",
            class = "bg-light",
            title = div(
              class = "mb-3 p-2 rounded",
              style = "background-color: white; border-left: 4px solid #3498DB;",
              div(
                class = "d-flex align-items-center",
                icon("bullseye", class = "text-info me-2"),
                strong("Foreground Selection Workflow")
              )
            ),

            # Accordion with organized sections
            bslib::accordion(
              id = "foreground_accordion",
              open = c("files", "config", "selection"), # Panels open by default

              # Step 1: Upload Files
              bslib::accordion_panel(
                title = div(
                  icon("upload", class = "me-2"),
                  "Step 1: Upload Files"
                ),
                value = "files",
                fileInput(
                  inputId = ns("geno_data"),
                  label = div(
                    icon("file-csv", class = "me-2 text-success"),
                    "Upload Genotypic Data (geno)"
                  ),
                  accept = c(".xlsx", ".xls", ".csv")
                ),
                fileInput(
                  inputId = ns("marker_info"),
                  label = div(
                    icon("file-csv", class = "me-2 text-info"),
                    "Upload Marker Info File"
                  ),
                  accept = c(".xlsx", ".xls", ".csv")
                ),
                textInput(
                  inputId = ns("sep"),
                  label = "Allele Separator",
                  value = ":"
                )
              ),

              # Step 2: Configure Columns
              bslib::accordion_panel(
                title = div(
                  icon("table-columns", class = "me-2"),
                  "Step 2: Configure Columns"
                ),
                value = "config",
                div(
                  class = "small text-muted mb-3",
                  icon("info-circle", class = "me-1"),
                  "Map your data columns to the required fields"
                ),
                selectInput(
                  inputId = ns("fore_marker_col"),
                  label = "Foreground Marker Column",
                  choices = NULL
                ),
                bslib::layout_columns(
                  col_widths = c(6, 6),
                  selectInput(
                    inputId = ns("fav_allele_col"),
                    label = "Favorable Allele:",
                    choices = NULL
                  ),
                  selectInput(
                    inputId = ns("alt_allele_col"),
                    label = "Alternative Allele:",
                    choices = NULL
                  )
                )
              ),

              # Step 3: Selection Parameters
              bslib::accordion_panel(
                title = div(
                  icon("sliders", class = "me-2"),
                  "Step 3: Selection Parameters"
                ),
                value = "selection",
                selectInput(
                  inputId = ns("select_type"),
                  label = "Selection Type",
                  choices = c("Homozygous" = "homo", "Heterozygous" = "hetero", "Both" = "both"),
                  selected = "homo"
                ),
                div(
                  class = "alert alert-info small mt-2",
                  icon("lightbulb", class = "me-1"),
                  "Choose the genotype selection criteria for your analysis"
                )
              )
            ),

            # Action button
            div(
              class = "mt-4 d-grid gap-2",
              actionButton(
                inputId = ns("run_analysis"),
                label = "Get Results",
                icon = icon("play", class = "me-2"),
                class = "btn-success btn-lg",
                style = "font-weight: 600; box-shadow: 0 4px 6px rgba(0,0,0,0.1);"
              )
            )
          ),

          # Main panel.
          bslib::input_switch(id = ns("config"), label = "Configure Plot", value = FALSE),
          conditionalPanel(
            condition = paste0('input["', ns("config"), '"] == true'),
            tagList(
              bslib::card(
                bslib::card_header(
                  div(
                    class = "d-flex align-items-center justify-content-between",
                    div(
                      icon("chart-bar", class = "me-2"),
                      strong("UpSet Plot Settings")
                    )
                  ),
                  class = "bg-primary text-white"
                ),
                bslib::card_body(
                  # Basic Settings Section
                  div(
                    class = "mb-4",
                    h5(class = "mb-3", icon("palette", class = "me-2"), "Basic Settings"),
                    bslib::layout_columns(
                      col_widths = c(4, 4, 4),
                      textInput(
                        inputId = ns("mainbar_y_label"),
                        label = "Main Bar Y-axis Label",
                        value = "Locus Intersection Size"
                      ),
                      textInput(
                        inputId = ns("sets_x_label"),
                        label = "Sets X-axis Label",
                        value = "Locus Size"
                      ),
                      numericInput(
                        inputId = ns("text_scale"),
                        label = "Text Scale",
                        value = 1.2,
                        min = 0.5,
                        max = 3,
                        step = 0.1
                      )
                    )
                  ),
                  hr(),

                  # Metadata Settings Section
                  div(
                    class = "mb-4",
                    h5(class = "mb-3", icon("tags", class = "me-2"), "Metadata Settings"),
                    bslib::layout_columns(
                      col_widths = c(4, 4, 4),
                      selectInput(
                        inputId = ns("plot_type"),
                        label = "Metadata Plot Type",
                        choices = c("text"),
                        selected = "text"
                      ),
                      numericInput(
                        inputId = ns("plot_assign"),
                        label = "Assign (Position)",
                        value = 8,
                        min = 1,
                        step = 1
                      ),
                      textInput(
                        inputId = ns("plot_column"),
                        label = "Metadata Column",
                        value = "locus"
                      )
                    ),
                    selectInput(
                      inputId = ns("plot_colors"),
                      label = "Colors",
                      selected = "darkblue",
                      choices = grDevices::colors()
                    )
                  )
                )
              )
            )
          ),
          # Results display section
          bslib::card(
            height = "auto",
            bslib::card_body(
              navset_card_tab(
                id = ns("results_tabs"),

                # Upset Plot Tab
                bslib::nav_panel(
                  title = "Upset Plot Visualization",
                  icon = icon("chart-line"),
                  # Upset Plot Output
                  plotOutput(ns("result_plot"), height = "800px"),
                  # Download plot card widget
                  bslib::card(card_footer(
                    fluidRow(
                      column(
                        3,
                        textInput(
                          inputId = ns("file_name"),
                          label = "Enter Filename",
                          value = "Upset_Plot"
                        )
                      ),
                      column(
                        3,
                        numericInput(
                          inputId = ns("width"),
                          label = "Set Plot Width",
                          value = 11, min = 1
                        )
                      ), column(
                        3,
                        numericInput(
                          inputId = ns("height"),
                          label = "Set Plot Height",
                          value = 10, min = 1
                        )
                      )
                    ),
                    downloadButton(
                      outputId = ns("download_plot"),
                      label = "Download Plot", class = "btn-success"
                    )
                  )),
                  hr(),

                  # Advanced Filtering Section
                  div(
                    h5(class = "mb-3", icon("filter", class = "me-2"), "Advanced Filtering"),
                    bslib::input_switch(
                      id = ns("configure"),
                      label = "Enable Line Filtering",
                      value = FALSE
                    ),

                    conditionalPanel(
                      condition = paste0('input["', ns("configure"), '"] == true'),
                      bslib::card(
                        class = "mt-3 bg-light",
                        bslib::card_body(
                          bslib::layout_columns(
                            col_widths = c(6, 6),

                            # Filter Controls
                            div(
                              h6(class = "mb-3", icon("gears", class = "me-2"), "Filter Criteria"),
                              selectInput(
                                inputId = ns("present"),
                                label = "SNPs to be Present",
                                choices = NULL,
                                multiple = TRUE,
                                width = '100%'
                              ),
                              selectInput(
                                inputId = ns("absent"),
                                label = "SNPs to be Absent",
                                choices = NULL,
                                multiple = TRUE,
                                width = '100%'
                              ),
                              div(
                                class = "d-grid gap-2 mt-3",
                                actionButton(
                                  inputId = ns("Search"),
                                  label = "Find Matching Lines",
                                  icon = icon("search"),
                                  class = "btn-primary btn-lg"
                                )
                              )
                            ),

                            # Display Results
                            div(
                              h6(class = "mb-3", icon("list", class = "me-2"), "Matching Genotypes"),
                              bslib::card(
                                style = "height: 280px; overflow-y: auto;",
                                bslib::card_body(
                                  verbatimTextOutput(outputId = ns("display_result"))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                ),

                # Data Table Tab
                bslib::nav_panel(
                  title = "Binary Matrix of Presence or Absence of favorable Alleles",
                  icon = icon("table"),
                  DT::DTOutput(ns("result_table"), height = "800px")
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
mod_ds_foreground_select_server <- function(id) {
  moduleServer(id, function(input, output, session) {
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
        selected = safe_grep_match(pattern = "marker", choices = cols)
      )

      updateSelectInput(session, "fav_allele_col",
        choices = cols,
        selected = safe_grep_match(pattern = "fav", choices = cols)
      )

      updateSelectInput(session,
        "alt_allele_col",
        choices = cols,
        selected = safe_grep_match(pattern = "alt", choices = cols)
      )
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
      # Update marker_info()
      marker_info_updated <- marker_info()[marker_info()[[input$fore_marker_col]] %in% colnames(geno_data()), ]

      result <- tryCatch({
        foreground_select(
          geno_data = geno_data(),
          fore_marker_info = marker_info_updated,
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
        shinyjs::delay(ms = 500, {
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

      tryCatch(
        {
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
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Error",
            text = paste("Upset plot failed:", e$message),
            type = "error"
          )
          return(NULL)
        }
      )
    })

    # Render upset plot
    output$result_plot <- renderPlot({
      req(upset_result())
      print(upset_result())
      # show success toast
      shinyWidgets::show_toast(
        title = "",
        text = "Upset plot has been successfully rendered!",
        type = "success",
        position = "bottom-end",
        timer = 3000
      )
    })

    # Download handler for the plot as pdf
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("upset_plot", input$file_name, ".pdf", sep = "")
      },
      content = function(file) {
        req(upset_result())
        grDevices::pdf(file, width = input$width, height = input$height) # Adjust dimensions as needed
        print(upset_result())
        grDevices::dev.off()
      }
    )

    # Update select input
    observe({
      req(result_data())

      updateSelectInput(session,
        inputId = "present",
        choices = colnames(result_data()),
        selected = character(0)
      )

      updateSelectInput(session,
        inputId = "absent",
        choices = colnames(result_data()),
        selected = character(0)
      )
    })



    found_lines <- eventReactive(input$Search, {
      # input$present , input$absent,
      req(result_data())

      find_lines(mat = result_data(), present = input$present, absent = input$absent)
    })



    # render result.
    output$display_result <- renderPrint({
      req(found_lines())
      paste0("line", found_lines())
    })
  })
}



## To be copied in the UI
# mod_ds_foreground_select_ui("ds_foreground_select_1")

## To be copied in the server
# mod_ds_foreground_select_server("ds_foreground_select_1")
