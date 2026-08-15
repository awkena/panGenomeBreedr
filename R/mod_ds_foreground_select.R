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
        title = "Foreground Selection",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            id = ns("sidebar"),
            width = 400,
            position = "left",
            class = "bg-light",
            title = "",

            # Accordion with organized sections
            bslib::accordion(
              id = "foreground_accordion",
              open = c("files", "config", "selection"), # Panels open by default

              ## Data Acquisition
              bslib::accordion_panel(
                title = div(
                  icon("upload", class = "me-2"),
                  tags$span(
                    "Data Acquisition",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "files",
                fileInput(
                  inputId = ns("geno_data"),
                  label = "Upload Marker Data",
                  accept = c(".xlsx", ".xls", ".csv")
                ),
                fileInput(
                  inputId = ns("marker_info"),
                  label = "Upload Marker Meta Data",
                  accept = c(".xlsx", ".xls", ".csv")
                ),
                textInput(
                  inputId = ns("sep"),
                  label = "Enter Allele Separator",
                  value = ":"
                )
              ),

              # Feature Mapping
              bslib::accordion_panel(
                title = div(
                  icon("project-diagram", class = "me-2"),
                  tags$span(
                    "Feature Mapping",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
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
                selectInput(
                  inputId = ns("sample_name"),
                  label = "Sample Name Column",
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
                  tags$span(
                    "Selection Parameters",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
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
          # div(
          #   id = ns("placeholder_ui"),
          #   class = "d-flex align-items-center justify-content-center text-muted h-100",
          #   style = "min-height: 600px; flex-direction: column;",
          #   icon("arrow-left", class = "mb-3", style = "font-size: 4rem; opacity: 0.5;"),
          #   h4("Upload and process your marker data to begin.", style = "opacity: 0.7;")
          # ),
          
          shinyjs::hidden(
            div(
              id = ns("analysis_ui"),
              bslib::navset_card_underline(
                id = ns("results_tabs"),
                
                bslib::nav_panel(
                  title = "Plot Configuration",
                  value = "config_tab",
                  icon = icon("sliders"),
                  bslib::card(
                    class = "shadow-sm border-0",
                    bslib::card_header(
                      div(
                        class = "d-flex align-items-center justify-content-between",
                        div(icon("chart-bar", class = "me-2"), strong("UpSet Plot Settings"))
                      ),
                      class = "bg-light"
                    ),
                    bslib::card_body(
                      div(
                        class = "mb-4",
                        h5(class = "mb-3", icon("palette", class = "me-2"), "Basic Settings"),
                        bslib::layout_columns(
                          col_widths = c(4, 4, 4),
                          textInput(inputId = ns("mainbar_y_label"), label = "Main Bar Y-axis Label", value = "Locus Intersection Size"),
                          textInput(inputId = ns("sets_x_label"), label = "Sets X-axis Label", value = "Locus Size"),
                          numericInput(inputId = ns("text_scale"), label = "Text Scale", value = 1.2, min = 0.5, max = 3, step = 0.1)
                        )
                      ),
                      hr(),
                      div(
                        class = "mb-4",
                        h5(class = "mb-3", icon("tags", class = "me-2"), "Metadata Settings"),
                        bslib::layout_columns(
                          col_widths = c(4, 4, 4),
                          selectInput(inputId = ns("plot_type"), label = "Metadata Plot Type", choices = c("text"), selected = "text"),
                          numericInput(inputId = ns("plot_assign"), label = "Assign (Position)", value = 8, min = 1, step = 1),
                          textInput(inputId = ns("plot_column"), label = "Metadata Column", value = "locus")
                        ),
                        selectInput(inputId = ns("plot_colors"), label = "Colors", selected = "darkblue", choices = grDevices::colors())
                      )
                    )
                  ),
                  div(
                    class = "d-flex justify-content-end mt-3",
                    actionButton(
                      inputId = ns("generate_plot"),
                      label = "Generate UpSet Plot",
                      icon = icon("play", class = "me-2"),
                      class = "btn-success btn-lg px-5",
                      style = "font-weight: 600;"
                    )
                  )
                ),

                bslib::nav_panel(
                  title = "UpSet Plot Viewer",
                  value = "plot_tab",
                  icon = icon("chart-area"),
                  div(
                    class = "mb-3 d-flex align-items-end gap-3",
                    div(
                      style = "flex: 1; max-width: 250px;",
                      textInput(inputId = ns("file_name"), label = "Enter Filename", value = "Upset_Plot", width = "100%")
                    ),
                    div(
                      style = "flex: 1; max-width: 150px;",
                      numericInput(inputId = ns("width"), label = "Plot Width", value = 11, min = 1, width = "100%")
                    ),
                    div(
                      style = "flex: 1; max-width: 150px;",
                      numericInput(inputId = ns("height"), label = "Plot Height", value = 10, min = 1, width = "100%")
                    ),
                    div(
                      class = "mb-3",
                      downloadButton(outputId = ns("download_plot"), label = "Export PDF", class = "btn-success")
                    )
                  ),
                  bslib::card(
                    class = "shadow-sm border-0",
                    full_screen = TRUE,
                    bslib::card_body(
                      shinycssloaders::withSpinner(
                        plotOutput(ns("result_plot"), height = "800px"),
                        type = 4, color = "#0dc5c1"
                      )
                    )
                  )
                ),
                
                bslib::nav_panel(
                  title = "Line Filtering",
                  value = "filter_tab",
                  icon = icon("filter"),
                  bslib::card(
                    class = "border-0 shadow-sm bg-light",
                    bslib::card_header(tags$h5(icon("filter", class = "me-2"), "Advanced Line Filtering", class = "mb-0")),
                    bslib::card_body(
                      bslib::layout_column_wrap(
                        width = 1/2,
                        gap = "20px",
                        bslib::card(
                          class = "border-0 shadow-sm",
                          bslib::card_header(tags$h6(icon("sliders", class = "me-2"), "Filter Criteria", class = "fw-bold mb-0")),
                          bslib::card_body(
                            selectInput(inputId = ns("present"), label = "SNPs to be Present", choices = NULL, multiple = TRUE, width = '100%'),
                            selectInput(inputId = ns("absent"), label = "SNPs to be Absent", choices = NULL, multiple = TRUE, width = '100%'),
                            div(
                              class = "d-grid gap-2 mt-3",
                              actionButton(inputId = ns("Search"), label = "Find Matching Lines", icon = icon("search"), class = "btn-primary btn-lg shadow-sm")
                            )
                          )
                        ),
                        bslib::card(
                          class = "border-0 shadow-sm",
                          bslib::card_header(
                            class = "d-flex justify-content-between align-items-center",
                            tags$h6(icon("list-check", class = "me-2"), "Matching Genotypes", class = "fw-bold mb-0"),
                            downloadButton(ns("download_results"), "Export CSV", class = "btn-sm btn-success border-0")
                          ),
                          bslib::card_body(
                            style = "height: 280px; overflow-y: auto; background-color: white; border: 1px solid #dee2e6; border-radius: 4px;",
                            verbatimTextOutput(outputId = ns("display_result"))
                          )
                        )
                      )
                    )
                  )
                )

                # bslib::nav_panel(
                #   title = "Binary Matrix",
                #   value = "table_tab",
                #   icon = icon("table"),
                #   bslib::card(
                #     class = "border-0 shadow-sm",
                #     bslib::card_body(
                #       DT::DTOutput(ns("result_table"), height = "800px")
                #     )
                #   )
                # )
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
    geno_data <- eventReactive(input$geno_data,{
      req(input$geno_data)
      read_mapfile(input$geno_data$datapath)
    })

    marker_info <- eventReactive(input$marker_info,{
      req(input$marker_info)
      read_mapfile(input$marker_info$datapath)
    })

    # Dynamically update select input choices when marker info is uploaded
    observe({
      cols <- colnames(marker_info()) # meta data colnames
      cols_2 <- colnames(geno_data()) # genotype data colnames

      updateSelectInput(session, "fore_marker_col",
        choices = cols,
        selected = safe_grep_match(pattern = "marker", choices = cols)
      )

      updateSelectizeInput(session, "sample_name",
                           server = T,
                           selected = cols_2[1],
                           choices = cols_2
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
        geno_data(), marker_info(),input$sample_name,
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

      ## Add rownames to geno data.
      df_geno <- geno_data()
      sample_col_data <- as.character(df_geno[[input$sample_name]])

      # Check for non-unique or missing names
      if (any(duplicated(sample_col_data)) || any(is.na(sample_col_data))) {

        # Make names unique
        repaired_names <- make.unique(ifelse(is.na(sample_col_data), "Unknown", sample_col_data))
        rownames(df_geno) <- repaired_names

      } else {

        rownames(df_geno) <- sample_col_data
      }

      result <- tryCatch({
        foreground_select(
          geno_data = df_geno,
          fore_marker_info = marker_info_updated,
          fore_marker_col = input$fore_marker_col,
          fav_allele_col = input$fav_allele_col,
          alt_allele_col = input$alt_allele_col,
          select_type = input$select_type,
          sep = check_sep(input$sep)
        )


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


    # Transition the UI once data is processed
    observeEvent(result_data(), {
      req(result_data())

      #shinyjs::hide("placeholder_ui")
      shinyjs::show("analysis_ui")
      bslib::toggle_sidebar(id = "sidebar", open = "closed")

      shinyWidgets::show_alert(
        title = "Analysis Complete",
        text = "Foreground selection matrix generated. You can now configure and view your UpSet plot.",
        type = "success"
      )
    })

    # # Render DT Table
    # output$result_table <- DT::renderDT({
    #   req(result_data())
    #   DT::datatable(result_data(), options = list(pageLength = 15))
    # })



     # Eagerly switch tabs when user clicks "Generate UpSet Plot"
    observeEvent(input$generate_plot, {
      req(result_data(), marker_info())
      bslib::nav_select("results_tabs", "plot_tab")
    })

    # Generating the upset plot.
    upset_result <- eventReactive(input$generate_plot, {
      req(
        result_data(),
        marker_info(),
        input$mainbar_y_label,
        input$sets_x_label,
        input$text_scale,
        input$plot_type,
        input$plot_assign,
        input$plot_column,
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



    # Empty reactive value container to hold results
    values <- reactiveValues(
      found_lines = NULL,
      df = NULL
    )

    observeEvent(input$Search, {
      req(result_data())

      # Find lines
      result <- find_lines(mat = result_data(), present = input$present, absent = input$absent)

      # Check if the result is NULL
      if (is.null(result) || length(result) == 0) {
        display_val <-  paste('\u2716 None Found')
        download_val <- character(0)
      } else {
        display_val <- result
        download_val <- result
      }

      # Store the result for display
      values$found_lines <- display_val

      values$df <- create_padded_df(
        snps_present = input$present,
        snps_absent = input$absent,
        sample_name = download_val
      ) # create a dataframe
    })


    # Print result.
    output$display_result <- renderPrint({
      req(values$found_lines)
       values$found_lines
    })



    # Download result.
    output$download_results <- downloadHandler(
      filename = function() {
        paste0("Find_lines_result", ".csv")
      },
      content = function(file) {
        req(values$df)
        # Write df as a csv file
        write.csv(values$df, file,row.names = F)
      })

  })
}



## To be copied in the UI
# mod_ds_foreground_select_ui("ds_foreground_select_1")

## To be copied in the server
# mod_ds_foreground_select_server("ds_foreground_select_1")
