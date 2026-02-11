#' Decision Support for Marker Assisted Backcrossing (UI)
#'
#' @description A Shiny module UI for analyzing marker data in backcross breeding programs.
#'
#' @param id Unique namespace ID
#'
#' @return A shiny::tagList containing the UI elements
#'
#' @importFrom shiny NS tagList fileInput radioButtons selectInput textInput
#' @importFrom shiny actionButton uiOutput icon plotOutput
#' @importFrom bslib navset_card_underline nav_panel layout_sidebar sidebar
#' @importFrom bslib card card_header card_body card_footer accordion
#' @importFrom bslib accordion_panel navset_card_tab input_switch
#' @noRd
#'
mod_ds_marker_ass_bac_ui <- function(id) {
  ns <- NS(id)
  tagList(
    navset_card_underline(
      nav_panel(
        title = "Decision Support for MABC",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            id = ns("sidebar"),
            width = 400,
            position = "left",
            class = "bg-light",

            # Accordion with organized sections
            bslib::accordion(
              id = "introgression_accordion",
              open = c("files", "mapfile", "settings", "qc"), # Panels open by default

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
                  inputId = ns("data_id"),
                  label = "Upload Genotype Data ",
                  placeholder = ".csv,  .xlsx",
                  multiple = FALSE,
                  accept = c(".csv", ".xlsx")
                ),
                layout_columns(
                  col_widths = c(6, 6),
                  selectInput(
                    inputId = ns("data_type"),
                    label = "Data Type:",
                    choices = c("Agriplex", "Kasp / DArTag"),
                    selected = "Agriplex",
                    width = "100%"
                  ),
                  numericInput(
                    inputId = ns("marker_start"),
                    label = "Marker Start:",
                    value = 7,
                    min = 1,
                    width = "100%"
                  )
                ),
                textInput(
                  inputId = ns("allele_sep"),
                  label = "Enter Allele Separator for Data Format",
                  value = " / ",
                  width = "100%"
                )
              ),

              ## Genomic Mapping Setup
              bslib::accordion_panel(
                title = div(
                  icon("map", class = "me-2"),
                  tags$span(
                    "Genomic Mapping Setup",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "mapfile",
                radioButtons(
                  inputId = ns("choice"),
                  label = "Do you have a map file?",
                  choices = c("Yes" = "yes", "Generate one" = "no"),
                  selected = "yes"
                ),

                # Dynamic rendering based on choice
                conditionalPanel(
                  condition = paste0('input["', ns("choice"), '"] == "yes"'),
                  tagList(
                    fileInput(
                      inputId = ns("mapfile"),
                      label = "Upload Map file",
                      accept = c(".csv", ".xlsx", ".xls"),
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("snp_id"),
                      label = "Select Column for SNP ID",
                      choices = NULL,
                      width = "100%"
                    )
                  )
                ),
                conditionalPanel(
                  condition = paste0('input["', ns("choice"), '"] == "no"'),
                  tagList(
                    textInput(
                      inputId = ns("prefix_marker"),
                      label = "Enter Prefix for Marker ID",
                      value = "S",
                      width = "100%"
                    ),
                    textInput(
                      inputId = ns("sep_marker"),
                      label = "Enter Separator for Marker ID",
                      value = "_",
                      width = "100%"
                    )
                  )
                )
                # div(
                #   class = "alert alert-info small mt-2",
                #   icon("lightbulb", class = "me-1"),
                #   "Map file links markers to their genomic positions"
                # )
              ),

              ## Selection Parameters
              bslib::accordion_panel(
                title = div(
                  icon("sliders", class = "me-2"),
                  tags$span(
                    "Selection Parameters",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "settings",
                selectizeInput(
                  inputId = ns("sample_id"),
                  label = "Select Sample ID Column",
                  choices = NULL,
                  width = "100%"
                ),
                selectizeInput(
                  inputId = ns("cluster_by"),
                  label = "Group By:",
                  choices = NULL,
                  width = "100%"
                ),
                uiOutput(outputId = ns("cluster_focus_ui")),
                bslib::layout_columns(
                  col_widths = c(6, 6),
                  div(
                    numericInput(
                      inputId = ns("rp"),
                      label = "Recurrent Parent Index:",
                      value = 1,
                      min = 1,
                      width = "100%"
                    ),
                    uiOutput(outputId = ns("recurrent_help"))
                  ),
                  div(
                    numericInput(
                      inputId = ns("dp"),
                      label = "Donor Parent Index:",
                      value = 3,
                      min = 1,
                      width = "100%"
                    ),
                    uiOutput(outputId = ns("donor_help"))
                  )
                )
              ),

              ## Genotypic Data Cleaning
              bslib::accordion_panel(
                title = div(
                  icon("filter", class = "me-2"),
                  tags$span(
                    "Genotypic Data Cleaning",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "qc",
                div(
                  class = "small text-muted mb-3",
                  icon("info-circle", class = "me-1"),
                  "Enable filters to improve data quality"
                ),
                bslib::input_switch(
                  id = ns("apply_par_poly"),
                  label = "Remove Monomorphic Parents",
                  value = TRUE
                ),
                bslib::input_switch(
                  id = ns("apply_par_miss"),
                  label = "Remove Missing Parent Data",
                  value = TRUE
                ),
                bslib::input_switch(
                  id = ns("apply_geno_good"),
                  label = "Apply Genotype Error Check",
                  value = TRUE
                ),
                bslib::input_switch(
                  id = ns("apply_par_homo"),
                  label = "Filter Heterozygous Parents",
                  value = TRUE
                )
              )
            ),

            # Action button to run analysis
            div(
              class = "mt-4 d-grid gap-2",
              actionButton(
                inputId = ns("config"),
                label = "Process Data",
                icon = icon("gears", class = "me-2"),
                class = "btn-primary btn-lg",
                style = "font-weight: 600;"
              )
            )
          ),

          # Main content area
          bslib::input_switch(id = ns("configure"), label = "Plot Configuration", value = T),
          conditionalPanel(
            condition = paste0('input["', ns("configure"), '"] == true'),
            card(
              fluidRow(
                # Card 1: Heatmap Configuration
                column(
                  width = 4,
                  bslib::card(
                    class = "shadow p",
                    max_height = "600px",
                    height = "100%",
                    bslib::card_header(tags$b("RPP Calculation Parameters"),
                      class = "bg-warning text-center",
                      style = "font-size:18px;"
                    ),
                    bslib::card_body(
                      selectInput(
                        inputId = ns("snp_ids"),
                        label = "Select Column for SNP ID",
                        choices = NULL,
                        width = "100%"
                      ),
                      selectInput(
                        inputId = ns("chr"),
                        label = "Select Column for Chromosome",
                        choices = NULL,
                        width = "100%"
                      ),
                      selectInput(
                        inputId = ns("chr_pos"),
                        label = "Select Column for Chromosome Position",
                        choices = NULL,
                        width = "100%"
                      ),
                      radioButtons(
                        inputId = ns("weight_rpp"),
                        label = "Use Weighted RPP?",
                        choices = c("Yes" = TRUE, "No" = FALSE),
                        selected = TRUE,
                        inline = TRUE
                      )
                    )
                  )
                ),

                # Card 2: BC Progenies RPP Plot Settings
                column(
                  width = 8,
                  bslib::card(
                    class = "shadow p",
                    max_height = "600px",
                    height = "100%",
                    bslib::card_header(tags$b("RPP Visualization Controls"),
                      class = "bg-info text-center",
                      style = "font-size:18px;"
                    ),
                    bslib::card_body(
                      fluidRow(
                        # Left Column - 4 widgets
                        column(
                          width = 6,
                          selectInput(
                            inputId = ns("rpp_col"),
                            label = "Select Total RPP Column",
                            choices = NULL,
                            width = "100%"
                          ),
                          selectInput(
                            inputId = ns("rpp_sample_id"),
                            label = "Select Sample ID Column",
                            choices = NULL,
                            width = "100%"
                          ),
                          numericInput(
                            inputId = ns("bc_gen"),
                            label = "BC Generation",
                            value = NULL,
                            min = 1,
                            width = "100%"
                          ),
                          numericInput(
                            inputId = ns("rpp_threshold"),
                            label = "Selection Threshold (%)",
                            value = 0.93,
                            min = 0,
                            max = 1,
                            width = "100%"
                          ),
                          selectInput(
                            inputId = ns("thresh_line_col"),
                            label = "Color of Threshold Line",
                            choices = grDevices::colors(),
                            selected = "firebrick",
                            width = "100%"
                          ),
                          bslib::input_switch(
                            id = ns("show_above_thresh"),
                            label = "Show Progenies Above Threshold",
                            value = FALSE
                          )
                        ),
                        # Right Column - 4 widgets
                        column(
                          width = 6,
                          selectInput(
                            inputId = ns("bar_col"),
                            label = "Set Bar Fill Color",
                            choices = grDevices::colors(),
                            selected = "cornflowerblue",
                            width = "100%"
                          ),
                          numericInput(
                            inputId = ns("alpha"),
                            label = "Adjust Bar Transparency",
                            value = 0.9,
                            min = 0,
                            max = 1,
                            step = 0.1,
                            width = "100%"
                          ),
                          numericInput(
                            inputId = ns("text_size"),
                            label = "Text Size",
                            value = 15,
                            min = 1,
                            width = "100%"
                          ),
                          numericInput(
                            inputId = ns("bar_width"),
                            label = "Set Bar Width",
                            value = 0.5,
                            min = 0.1,
                            width = "100%",
                            step = 0.1
                          ),
                          numericInput(
                            inputId = ns("aspect_ratio"),
                            label = "Set Aspect Ratio of Barplot",
                            value = 0.5,
                            min = 0.1,
                            width = "100%",
                            step = 0.1
                          ),
                          numericInput(
                            inputId = ns("text_scale_fct"),
                            label = "Set Text Size Scaling Factor",
                            value = 0.1,
                            min = 0.1,
                            width = "100%",
                            step = 0.1
                          )
                        )
                      )
                    )
                  )
                )
              ),
              actionButton(
                inputId = ns("generate_rpp_plot"),
                label = "Generate Plot",
                icon = icon("play", class = "me-2"),
                class = "btn-success btn-lg",
                style = "font-weight: 600;"
              )
            )
          ),

          # Results Section
          fluidRow(
            column(
              width = 12,
              bslib::accordion(
                bslib::accordion_panel(
                  title = "Recurrent Parent Percentage (%) Plot",
                  icon = icon("chart-line"),
                  div(
                    plotOutput(
                      outputId = ns("rpp_bar"),
                      width = "100%",
                      height = "600px"
                    ),
                    bslib::card(card_footer(
                      fluidRow(
                        column(
                          3,
                          textInput(
                            inputId = ns("file_name2"),
                            label = "Enter Filename",
                            value = "rpp_barplot"
                          )
                        ),
                        column(
                          3,
                          numericInput(
                            inputId = ns("width2"),
                            label = "Set Plot Width",
                            value = 8, min = 1
                          )
                        ), column(
                          3,
                          numericInput(
                            inputId = ns("height2"),
                            label = "Set Plot Height",
                            value = 6, min = 1
                          )
                        )
                      ),
                      downloadButton(
                        outputId = ns("download_plot2"),
                        label = "Download Plot", class = "btn-success"
                      )
                    ))
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

#' Decision Support for Marker Assisted Backcrossing (Server)
#'
#' @description Server logic for the MABC analysis module
#'
#' @param id Unique namespace ID
#'
#' @importFrom shiny moduleServer observeEvent req reactive renderUI observe
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @noRd
mod_ds_marker_ass_bac_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # # Dynamic ui based on choice of individual.
    observe({
      if (input$choice == "no") {
        showModal(
          modalDialog(
            title = div(
              style = "display: flex; align-items: center;",
              icon("circle-info", class = "text-primary me-2"),
              tags$h4("Marker Formatting Guide", class = "mb-0", style = "font-weight: 700;")
            ),
            tagList(
              p("To ensure your markers are parsed correctly into the map file, please follow this naming convention:",
                class = "text-muted mb-4"
              ),

              # Using a list with custom spacing
              tags$ul(
                class = "list-group list-group-flush mb-4",
                tags$li(
                  class = "list-group-item border-0 ps-0",
                  icon("chevron-right", class = "text-primary me-2 small"),
                  "A common prefix (e.g., ", tags$code("S"), ")"
                ),
                tags$li(
                  class = "list-group-item border-0 ps-0",
                  icon("chevron-right", class = "text-primary me-2 small"),
                  "Chromosome number (e.g., ", tags$code("1"), ")"
                ),
                tags$li(
                  class = "list-group-item border-0 ps-0",
                  icon("chevron-right", class = "text-primary me-2 small"),
                  "A separator character (e.g., ", tags$code("_"), ")"
                ),
                tags$li(
                  class = "list-group-item border-0 ps-0",
                  icon("chevron-right", class = "text-primary me-2 small"),
                  "Position number (e.g., ", tags$code("101"), ")"
                )
              ),

              # Example to user, highlighted in a Card-style div
              div(
                class = "p-3 rounded-3 bg-light border",
                style = "border-left: 5px solid #0d6efd !important;",
                tags$b("Example: ", class = "text-primary"),
                tags$code("S1_101"),
                p(
                  class = "small text-muted mt-2 mb-0",
                  "This indicates Prefix (S), Chromosome (1), and Position (101)."
                )
              )
            ),
            easyClose = FALSE,
            size = "m",
            footer = div(
              class = "d-flex justify-content-end",
              actionButton(
                inputId = ns("close_modal"),
                label = "Understood",
                class = "btn-primary px-4 shadow-sm",
                style = "border-radius: 8px;",
                onclick = "Shiny.setInputValue('close_modal', true, {priority: 'event'})"
              )
            )
          )
        )
      }
    })

    # Close info modal
    observeEvent(input$close_modal, {
      removeModal()
    })

    # Read and validate uploaded CSV file
    validated_data <- eventReactive(input$data_id, {
      req(input$data_id)

      read_mapfile(filepath = input$data_id$datapath)
    })


    # Get column names from the validated data
    data_colnames <- reactive({
      req(validated_data())
      colnames(validated_data())
    })


    # Populate Batch and Genotype column selectors dynamically
    observe({
      req(data_colnames())

      # Update batch column selector
      updateSelectizeInput(
        session,
        inputId = "cluster_by",
        choices = c("None", data_colnames()), server = T
      )

      # Update genotype column selector
      updateSelectizeInput(
        session,
        inputId = "sample_id",
        choices = data_colnames(),
        selected = safe_grep_match(
          pattern = "genotype",
          choices = data_colnames()
        ), server = T
      )
    })


    # cluster logic
    output$cluster_focus_ui <- renderUI({
      # Only show the second dropdown if a factor is selected
      req(input$cluster_by, validated_data())
      if (input$cluster_by == "None") {
        return(NULL)
      }

      # Get the unique values
      choices_to_focus <- unique(validated_data()[[input$cluster_by]])

      selectInput(
        inputId = ns("cluster_focus"),
        label = paste("Select", input$cluster_by, "to focus on:"),
        choices = choices_to_focus
      )
    })


    # Update the validated data based on the cluster selected.
    sub_validated_data <- reactive({
      req(validated_data(), input$cluster_by)

      if (!input$cluster_by %in% colnames(validated_data()) || input$cluster_by == "None") {
        return(validated_data())
      } else {
        req(input$cluster_focus)
        return(validated_data()[validated_data()[[input$cluster_by]] == input$cluster_focus, ])
      }
    })


    ## Help text for donor and reccurent parents
    # Display the Donor Name
    output$donor_help <- renderUI({
      if (input$cluster_by != "None") {
        req(input$dp, input$sample_id, sub_validated_data())
        name <- sub_validated_data()[[input$sample_id]][input$dp]
        helpText(strong("Selected: "), span(name, style = "color: #e67e22;"))
      } else {
        # if user's data cannot be partitioned into groups
        req(validated_data(), input$sample_id, input$dp)
        name <- validated_data()[input$dp, input$sample_id]
        helpText(strong("Selected: "), span(name, style = "color: #0275d8; font-size: 1.1em;"))
      }
    })

    #  Display the Recurrent Name
    output$recurrent_help <- renderUI({
      if (input$cluster_by != "None") {
        req(input$rp, input$sample_id, sub_validated_data())
        name <- sub_validated_data()[[input$sample_id]][input$rp]
        helpText(strong("Selected: "), span(name, style = "color: #0275d8; font-size: 1.1em;"))
      } else {
        # if user's data cannot be partitioned into groups
        req(validated_data(), input$sample_id, input$rp)
        name <- validated_data()[input$rp, input$sample_id]
        helpText(strong("Selected: "), span(name, style = "color: #0275d8; font-size: 1.1em;"))
      }
    })


    # Read map file if user has.
    map_file <- reactive({
      req(input$mapfile)
      read_mapfile(filepath = input$mapfile$datapath)
    })

    # Populate field for snp id column based on mapfile
    observe({
      req(map_file())
      updateSelectInput(session,
        inputId = "snp_id",
        choices = colnames(map_file()),
        selected = safe_grep_match(
          pattern = "id",
          choices = colnames(map_file())
        )
      )
    })


    # process data.
    Result <- eventReactive(input$config, {
      req(
        sub_validated_data(), input$marker_start, input$sample_id,
        input$allele_sep, input$rp, input$dp, input$choice
      )

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Processing... Please wait."
      )

      result <- tryCatch({
        proc_nd_map_func(
          data = sub_validated_data(),
          marker_start = input$marker_start,
          sample_id = input$sample_id,
          marker_sep = if (input$choice == "no") input$sep_marker else NULL,
          apply_par_poly = input$apply_par_poly,
          apply_par_miss = input$apply_par_miss,
          apply_geno_good = input$apply_geno_good,
          apply_par_homo = input$apply_par_homo,
          snp_id = if (input$choice == "yes") input$snp_id else NULL,
          calls_sep = check_sep(input$allele_sep),
          data_type = if (input$data_type == "Agriplex") "agriplex" else if (input$data_type == "Kasp / DArTag") "kasp",
          rp = input$rp,
          dp = input$dp,
          Prefix = if (input$choice == "no") input$prefix_marker else NULL,
          feedback = input$choice,
          mapfile_path = if (!is.null(map_file())) map_file() else NULL
        )
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Error",
          text = paste("An error occurred while processing the data:", e$message),
          type = "error"
        )
        return(NULL)
      }, finally = {
        shinyjs::delay(ms = 1000, {
          shinybusy::remove_modal_spinner()
        })
      })

      result
    })

    # Update user on next step
    observeEvent(Result(), {
      req(Result()) # Only proceed if Result is not NULL

      # bslib::sidebar_toggle(id = 'sidebar' , open = "closed")
      shinyWidgets::show_alert(
        title = "Data Processed!",
        text = "Your genotype data is ready. Please proceed to 'Plot Configuration'  to generate your heatmap.",
        type = "success",
        btn_labels = "Got it!",
        closeOnClickOutside = FALSE,
        showCloseButton = TRUE
      )
    })


    # Observer for updating map-related inputs
    observe({
      req(Result())

      # Update SNP ID selection
      updateSelectInput(session,
        inputId = "snp_ids",
        choices = colnames(Result()$mapfile),
        selected = safe_grep_match(
          pattern = "id",
          choices = colnames(Result()$mapfile)
        )
      )

      # Update Chromosome selection
      updateSelectInput(session,
        inputId = "chr",
        choices = colnames(Result()$mapfile),
        selected = safe_grep_match(
          pattern = "chr",
          choices = colnames(Result()$mapfile)
        )
      )

      # Update Position selection
      updateSelectInput(session,
        inputId = "chr_pos",
        choices = colnames(Result()$mapfile),
        selected = safe_grep_match(
          pattern = "pos",
          choices = colnames(Result()$mapfile)
        )
      )
    })


    calc_rpp_bc_result <- reactive({
      req(
        Result(), input$chr_pos, input$chr, input$weight_rpp,
        input$snp_ids, sub_validated_data(), input$rp, input$sample_id
      )

      # Get recurrent parent name using index index
      parent_name <- sub_validated_data()[input$rp, input$sample_id]

      # Ensure the extracted name is present in the result's rowname
      req(parent_name %in% rownames(Result()$proc_kasp_f))

      calc_rpp_bc(
        x = Result()$proc_kasp_f,
        map_file = Result()$mapfile,
        map_chr = input$chr,
        map_pos = input$chr_pos,
        map_snp_ids = input$snp_ids,
        rp = parent_name,
        na_code = -5,
        weighted = input$weight_rpp
      )
    })

    # Update UI to enable plot generation
    observe({
      data <- calc_rpp_bc_result()
      req(data)

      cols <- colnames(data)

      updateSelectInput(session, "rpp_col",
        choices = cols,
        selected = safe_grep_match("total_rpp", cols)
      )

      updateSelectInput(session, "rpp_sample_id",
        choices = cols,
        selected = safe_grep_match("sample_id", cols)
      )
    })

    #  Generate the Barplot when user clicks the button
    rpp_barplot_result <- eventReactive(input$generate_rpp_plot, {
      req(calc_rpp_bc_result(), input$rpp_col, input$rpp_sample_id)

      tryCatch(
        {
          rpp_barplot(
            rpp_df = calc_rpp_bc_result(),
            rpp_sample_id = input$rpp_sample_id,
            rpp_col = input$rpp_col,
            rpp_threshold = input$rpp_threshold,
            text_size = input$text_size,
            text_scale_fct = input$text_scale_fct,
            alpha = input$alpha,
            bar_width = input$bar_width,
            aspect_ratio = input$aspect_ratio,
            bar_col = input$bar_col,
            thresh_line_col = input$thresh_line_col,
            show_above_thresh = input$show_above_thresh,
            bc_gen = if (is.null(input$rpp_threshold)) input$bc_gen else NULL,
            pdf = FALSE
          )
        },
        error = function(e) {
          showNotification(paste("Plot error:", e$message), type = "error")
          return(NULL)
        }
      )
    })


    # plot it.
    output$rpp_bar <- renderPlot({
      tryCatch(
        {
          plot_obj <- rpp_barplot_result()
          req(plot_obj)

          if (inherits(plot_obj, "ggplot")) {
            print(plot_obj)

            shinyjs::delay(ms = 1500, {
              shinyWidgets::show_toast(
                title = "",
                text = "Plot rendered successfully",
                type = "success"
              )
            })
          } else {
            plot_obj
          }
        },
        error = function(e) {
          return(NULL)
          shinyWidgets::show_alert(
            title = "Error",
            text = paste("An error occurred while generating plot", e$message),
            type = "error"
          )
        }
      )
    })


    # Download rpp plot
    output$download_plot2 <- downloadHandler(
      filename = function() {
        req(input$file_name2)
        paste0(input$file_name2, ".pdf")
      },
      content = function(file) {
        req(rpp_barplot_result())

        tryCatch(
          {
            grDevices::pdf(file,
              width = input$width2,
              height = input$height2,
              onefile = FALSE
            )

            # Just print the single plot
            print(rpp_barplot_result())

            grDevices::dev.off()
          },
          error = function(e) {
            shinyWidgets::show_toast(
              title = "Error",
              type = "error",
              text = paste("Download failed:", e$message),
              timer = 5000
            )
          }
        )
      }
    )
  })
}

## To be copied in the UI
# mod_ds_marker_ass_bac_ui("ds_marker_ass_bac_1")

## To be copied in the server
# mod_ds_marker_ass_bac_server("ds_marker_ass_bac_1")
