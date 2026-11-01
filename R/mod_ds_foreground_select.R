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

  # Accordion section title with an icon
  section_title <- function(icon_name, text) {
    div(icon(icon_name, class = "me-2"), tags$span(text, class = "fw-semibold"))
  }

  tagList(
    navset_card_underline(
      nav_panel(
        title = "Foreground Selection",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            id = ns("sidebar"),
            width = 350,
            position = "left",
            class = "bg-light",

            # SIDEBAR -------------------------------------------------------

            bslib::accordion(
              id = "foreground_accordion",
              open = "files", # Only the first panel is open by default

              # Upload files
              bslib::accordion_panel(
                title = section_title("upload", "Data Acquisition"),
                value = "files",
                fileInput(
                  inputId = ns("geno_data"),
                  label = label_help(
                    "Genotype Calls",
                    "One row per line and one column per marker, with calls such as A:G."
                  ),
                  placeholder = ".csv, .xlsx",
                  accept = c(".xlsx", ".xls", ".csv")
                ),
                fileInput(
                  inputId = ns("marker_info"),
                  label = label_help(
                    "Marker Information",
                    "One row per trait-predictive marker, giving its favourable and alternative allele."
                  ),
                  placeholder = ".csv, .xlsx",
                  accept = c(".xlsx", ".xls", ".csv")
                ),
                example_data_link(ns("load_example")),
                uiOutput(ns("example_badge")),
                tags$details(
                  class = "small text-muted mb-3",
                  tags$summary("Data format"),
                  tags$ul(
                    class = "ps-3 mt-2 mb-0",
                    tags$li("Genotype calls: a column of line names, then one column per marker (e.g. SNP1 = A:G)."),
                    tags$li("Marker information: marker name (matching the genotype column names), favourable allele, alternative allele, and optionally a locus name used to label the plot."),
                    tags$li("Markers in the information file that are not in the genotype file are skipped.")
                  )
                ),
                textInput(
                  inputId = ns("sep"),
                  label = label_help("Allele Separator", "The character between the two alleles of a call, e.g. ':' in A:G."),
                  value = ":"
                )
              ),

              # Map columns
              bslib::accordion_panel(
                title = section_title("project-diagram", "Map Columns"),
                value = "config",
                selectInput(
                  inputId = ns("sample_name"),
                  label = "Line Name Column (genotype file)",
                  choices = NULL
                ),
                selectInput(
                  inputId = ns("fore_marker_col"),
                  label = "Marker Column (information file)",
                  choices = NULL
                ),
                bslib::layout_columns(
                  col_widths = c(6, 6),
                  selectInput(
                    inputId = ns("fav_allele_col"),
                    label = "Favourable Allele",
                    choices = NULL
                  ),
                  selectInput(
                    inputId = ns("alt_allele_col"),
                    label = "Alternative Allele",
                    choices = NULL
                  )
                )
              ),

              # Homo, het or either
              bslib::accordion_panel(
                title = section_title("sliders", "Selection"),
                value = "selection",
                radioButtons(
                  inputId = ns("select_type"),
                  label = label_help(
                    "A line carries the favourable allele when it is",
                    "Homozygous: two copies (e.g. A:A), the allele is fixed. Heterozygous: one copy (e.g. A:G), it can still segregate. Both: either one."
                  ),
                  choices = c("Homozygous" = "homo", "Heterozygous" = "hetero", "Either" = "both"),
                  selected = "homo"
                )
              )
            ),

            run_button(ns("run_analysis"), "Process Data", icon_name = "gears")
          ),

          # PLACEHOLDER ---------------------------------------------------

          # Shown until the data has been processed
          div(
            id = ns("placeholder_ui"),
            class = "text-center text-muted py-5",
            div(class = "h1 mb-3", icon("bullseye")),
            tags$p(
              "Upload genotype calls and marker information in the sidebar, or load the example data, then click",
              tags$b("Process Data"), "."
            ),
            tags$p(
              class = "small",
              "You'll then see which lines carry the favourable allele at each target locus, and which lines combine the loci you want."
            )
          ),

          # RESULTS -------------------------------------------------------

          # Hidden until the data has been processed
          shinyjs::hidden(
            div(
              id = ns("analysis_ui"),
              bslib::navset_card_underline(
                id = ns("results_tabs"),

                # Selection results tab
                bslib::nav_panel(
                  title = "Selection Results",
                  value = "table_tab",
                  icon = icon("table"),
                  ui_result_card(
                    "Favourable alleles per line",
                    full_screen = TRUE,
                    actions = downloadButton(ns("download_matrix"), "CSV", class = "btn-sm btn-success"),
                    uiOutput(ns("matrix_summary")),
                    DT::DTOutput(ns("result_table"))
                  )
                ),

                # UpSet plot tab
                bslib::nav_panel(
                  title = "UpSet Plot",
                  value = "plot_tab",
                  icon = icon("chart-column"),
                  ui_result_card(
                    "Lines sharing each combination of favourable alleles",
                    full_screen = TRUE,
                    actions = tagList(
                      plot_options_button(
                        textInput(inputId = ns("mainbar_y_label"), label = "Bar Axis Label", value = "Locus Intersection Size"),
                        textInput(inputId = ns("sets_x_label"), label = "Set Size Axis Label", value = "Locus Size"),
                        bslib::layout_columns(
                          numericInput(inputId = ns("text_scale"), label = "Text Scale", value = 1.2, min = 0.5, max = 3, step = 0.1),
                          numericInput(inputId = ns("plot_assign"), label = "Locus Label Width", value = 8, min = 1, step = 1)
                        ),
                        colourpicker::colourInput(inputId = ns("plot_colors"), label = "Locus Label Colour", value = "darkblue", returnName = TRUE),
                        actionButton(
                          inputId = ns("redraw_plot"),
                          label = "Apply",
                          icon = icon("rotate"),
                          class = "btn-success w-100"
                        ),
                        title = "Plot appearance"
                      ),
                      plot_download_button(
                        ns, "download_plot",
                        file_name = "upset_plot",
                        formats = c("PDF" = "pdf"),
                        width = 11, height = 10
                      )
                    ),
                    tags$p(
                      class = "small text-muted",
                      "Each bar counts the lines that carry exactly the favourable alleles marked by the dots beneath it. The bars on the left count the lines carrying each locus."
                    ),
                    shinycssloaders::withSpinner(
                      plotOutput(ns("result_plot"), height = "700px"),
                      type = 4, color = "#0dc5c1"
                    )
                  )
                ),

                # Find lines tab
                bslib::nav_panel(
                  title = "Find Lines",
                  value = "filter_tab",
                  icon = icon("filter"),
                  bslib::layout_columns(
                    col_widths = bslib::breakpoints(sm = 12, lg = c(5, 7)),
                    ui_result_card(
                      "Combination you want",
                      selectizeInput(
                        inputId = ns("present"),
                        label = label_help("Favourable allele present at", "Lines must carry the favourable allele at every marker listed here."),
                        choices = NULL, multiple = TRUE, width = "100%"
                      ),
                      selectizeInput(
                        inputId = ns("absent"),
                        label = label_help("Favourable allele absent at", "Lines must not carry the favourable allele at any marker listed here."),
                        choices = NULL, multiple = TRUE, width = "100%"
                      ),
                      run_button(ns("Search"), "Find Matching Lines", icon_name = "search")
                    ),
                    ui_result_card(
                      "Matching lines",
                      actions = downloadButton(ns("download_results"), "CSV", class = "btn-sm btn-success"),
                      uiOutput(ns("display_result"))
                    )
                  )
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

    # DATA UPLOAD -----------------------------------------------------------

    # Uploaded (or example) data
    geno_data <- reactiveVal(NULL)
    marker_info <- reactiveVal(NULL)
    example_active <- reactiveVal(FALSE)

    # Genotype file
    observeEvent(input$geno_data, {
      geno_data(read_mapfile(input$geno_data$datapath))
      example_active(FALSE)
    })

    # Marker information file
    observeEvent(input$marker_info, {
      marker_info(read_mapfile(input$marker_info$datapath))
      example_active(FALSE)
    })

    # EXAMPLE DATA ----------------------------------------------------------

    # Vignette example: 4 lines, 3 trait-predictive markers
    observeEvent(input$load_example, {
      geno_data(data.frame(
        Line = paste0("Line", 1:4),
        SNP1 = c("A:A", "A:G", "G:G", "A:A"),
        SNP2 = c("C:C", "C:T", "T:T", "C:T"),
        SNP3 = c("G:G", "G:G", "A:G", "A:A")
      ))
      marker_info(data.frame(
        qtl_markers = paste0("SNP", 1:3),
        locus_name = paste0("loc", 1:3),
        fav_alleles = c("A", "C", "G"),
        alt_alleles = c("G", "T", "A")
      ))
      updateTextInput(session, "sep", value = ":")
      example_active(TRUE)
    })

    # Clear the example
    observeEvent(input$clear_example, {
      geno_data(NULL)
      marker_info(NULL)
      example_active(FALSE)
    })

    # "Using example data" badge
    output$example_badge <- renderUI({
      req(example_active())
      example_data_badge(
        "vignette foreground example (4 lines, 3 markers)",
        clear_id = ns("clear_example")
      )
    })

    # COLUMN DROPDOWNS ------------------------------------------------------

    # Fill the dropdowns once both files are in
    observe({
      req(marker_info(), geno_data())
      cols <- colnames(marker_info())
      cols_2 <- colnames(geno_data())

      updateSelectInput(session, "fore_marker_col",
        choices = cols,
        selected = safe_grep_match(pattern = "marker", choices = cols)[1]
      )

      updateSelectizeInput(session, "sample_name",
                           server = TRUE,
                           selected = cols_2[1],
                           choices = cols_2
      )

      updateSelectInput(session, "fav_allele_col",
        choices = cols,
        selected = safe_grep_match(pattern = "fav", choices = cols)[1]
      )

      updateSelectInput(session,
        "alt_allele_col",
        choices = cols,
        selected = safe_grep_match(pattern = "alt", choices = cols)[1]
      )
    })

    # RUN SELECTION ---------------------------------------------------------

    # Run foreground_select() on Process Data
    result_data <- eventReactive(input$run_analysis, {
      if (!validate_required_inputs(list(
        "Genotype Calls" = geno_data(),
        "Marker Information" = marker_info(),
        "Line Name Column" = input$sample_name,
        "Marker Column" = input$fore_marker_col,
        "Favourable Allele" = input$fav_allele_col,
        "Alternative Allele" = input$alt_allele_col,
        "Allele Separator" = input$sep
      ))) {
        return(NULL)
      }

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Running foreground analysis... Please wait."
      )
      # Keep only the markers that are in the genotype file
      marker_info_updated <- marker_info()[marker_info()[[input$fore_marker_col]] %in% colnames(geno_data()), ]

      # Line names become rownames
      df_geno <- geno_data()
      sample_col_data <- as.character(df_geno[[input$sample_name]])

      # Make missing or repeated names unique
      if (any(duplicated(sample_col_data)) || any(is.na(sample_col_data))) {
        repaired_names <- make.unique(ifelse(is.na(sample_col_data), "Unknown", sample_col_data))
        rownames(df_geno) <- repaired_names
      } else {
        rownames(df_geno) <- sample_col_data
      }

      result <- tryCatch({
        if (nrow(marker_info_updated) == 0) {
          stop("none of the markers in the marker information file are columns of the genotype file.")
        }
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
          title = "Couldn't run foreground selection",
          text = paste0(e$message, " Check the column choices under Map Columns and the allele separator."),
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

    # Swap placeholder for results once data is processed
    observeEvent(result_data(), {
      req(result_data())

      shinyjs::hide("placeholder_ui")
      shinyjs::show("analysis_ui")
      bslib::nav_select("results_tabs", "table_tab")
      bslib::toggle_sidebar(id = "sidebar", open = "closed")

      shinyWidgets::show_alert(
        title = "Selection complete",
        text = "See which lines carry the favourable alleles, the UpSet plot of their combinations, and Find Lines to pick a combination.",
        type = "success"
      )
    })

    # SELECTION RESULTS -----------------------------------------------------

    # Selection type in words, for the summary line
    select_label <- reactive({
      switch(input$select_type,
        homo = "homozygous for",
        hetero = "heterozygous for",
        both = "carrying at least one copy of"
      )
    })

    # Summary line above the table
    output$matrix_summary <- renderUI({
      mat <- result_data()
      req(mat)
      n_fav <- rowSums(mat == 1, na.rm = TRUE)
      k <- ncol(mat)
      tags$p(
        class = "small text-muted",
        tags$b(sum(n_fav == k)), " of ", nrow(mat), " lines are ", select_label(),
        " the favourable allele at all ", k, if (k == 1) " marker" else " markers",
        "; ", tags$b(sum(n_fav == 0)), " carry it at none.",
        HTML(" &#10003; = carried, &ndash; = not carried, ? = missing call.")
      )
    })

    # Results table: one row per line, a tick per locus carried
    output$result_table <- DT::renderDT({
      mat <- result_data()
      req(mat)
      n_fav <- rowSums(mat == 1, na.rm = TRUE)
      sym_yes <- '<span style="color:#1e7e34;font-weight:bold;">&#10003;</span>'
      sym_no <- '<span style="color:#adb5bd;">&ndash;</span>'
      sym_missing <- '<span style="color:#c0392b;font-weight:bold;">?</span>'
      shown <- as.data.frame(lapply(mat, function(x) {
        ifelse(is.na(x), sym_missing, ifelse(x == 1, sym_yes, sym_no))
      }), check.names = FALSE)
      shown <- data.frame(
        line = rownames(mat),
        loci_carried = paste0(n_fav, " of ", ncol(mat)),
        shown,
        check.names = FALSE
      )
      shown <- shown[order(-n_fav), , drop = FALSE]
      DT::datatable(
        shown,
        rownames = FALSE,
        escape = FALSE,
        colnames = c("Line", "Loci carried", colnames(mat)),
        options = list(
          pageLength = 15, scrollX = TRUE,
          columnDefs = list(list(className = "dt-center", targets = seq_len(ncol(mat)) + 1))
        )
      )
    })

    # Results table as csv
    output$download_matrix <- downloadHandler(
      filename = function() "foreground_selection.csv",
      content = function(file) {
        mat <- result_data()
        req(mat)
        out <- data.frame(line = rownames(mat), loci_carried = rowSums(mat == 1, na.rm = TRUE),
                          mat, check.names = FALSE)
        utils::write.csv(out, file, row.names = FALSE)
      }
    )

    # UPSET PLOT ------------------------------------------------------------

    # Drawn after processing, redrawn by Apply in Appearance
    upset_result <- eventReactive(list(result_data(), input$redraw_plot), {
      req(result_data(), marker_info())

      tryCatch(
        {
          run_upset_plot(
            foreground_matrix = result_data(),
            marker_info = marker_info(),
            mainbar_y_label = input$mainbar_y_label,
            sets_x_label = input$sets_x_label,
            text_scale = input$text_scale,
            plot_type = "text",
            assign = input$plot_assign,
            column = "locus",
            colors = input$plot_colors
          )
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Couldn't draw the UpSet plot",
            text = e$message,
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

    # Plot as pdf
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(input$download_plot_name, ".pdf")
      },
      content = function(file) {
        req(upset_result())
        grDevices::pdf(file, width = input$download_plot_width, height = input$download_plot_height, onefile = FALSE)
        print(upset_result())
        grDevices::dev.off()
      }
    )

    # FIND LINES ------------------------------------------------------------

    # Fill the present / absent lists with the markers
    observe({
      req(result_data())

      updateSelectizeInput(session,
        inputId = "present",
        choices = colnames(result_data()),
        selected = character(0)
      )

      updateSelectizeInput(session,
        inputId = "absent",
        choices = colnames(result_data()),
        selected = character(0)
      )
    })

    # Search results
    values <- reactiveValues(
      found_lines = NULL,
      df = NULL
    )

    # Find Matching Lines
    observeEvent(input$Search, {
      req(result_data())

      if (length(intersect(input$present, input$absent))) {
        shinyWidgets::show_toast(
          title = "",
          text = "A marker can't be in both lists. Remove it from one of them.",
          type = "warning"
        )
        return()
      }

      result <- find_lines(mat = result_data(), present = input$present, absent = input$absent)

      values$found_lines <- if (is.null(result)) character(0) else result

      values$df <- create_padded_df(
        snps_present = input$present,
        snps_absent = input$absent,
        sample_name = values$found_lines
      )
    })

    # Matching lines as chips
    output$display_result <- renderUI({
      if (is.null(values$found_lines)) {
        return(tags$p(
          class = "small text-muted",
          "Choose the markers, then click Find Matching Lines. Leave both boxes empty to list every line."
        ))
      }
      found <- values$found_lines
      tagList(
        tags$p(
          class = "small text-muted",
          tags$b(length(found)), " of ", nrow(result_data()), " lines match."
        ),
        if (length(found)) {
          # Neutral chips so results don't look like buttons; scroll when the list is long
          div(
            class = "d-flex flex-wrap gap-2",
            style = "max-height: 480px; overflow-y: auto;",
            lapply(found, function(x) {
              tags$span(
                class = "badge border fw-normal",
                style = "font-size: 0.85rem; background-color: #f1f3f5; color: #212529;",
                x
              )
            })
          )
        } else {
          tags$p(icon("circle-xmark", class = "text-danger me-1"), "No line has this combination.")
        }
      )
    })

    # Matching lines as csv
    output$download_results <- downloadHandler(
      filename = function() {
        paste0("Find_lines_result", ".csv")
      },
      content = function(file) {
        req(values$df)
        utils::write.csv(values$df, file, row.names = FALSE)
      })
  })
}

## To be copied in the UI
# mod_ds_foreground_select_ui("ds_foreground_select_1")

## To be copied in the server
# mod_ds_foreground_select_server("ds_foreground_select_1")
