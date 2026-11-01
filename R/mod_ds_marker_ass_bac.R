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

  # Accordion section title with an icon
  section_title <- function(icon_name, text) {
    div(icon(icon_name, class = "me-2"), tags$span(text, class = "fw-semibold"))
  }

  tagList(
    navset_card_underline(
      nav_panel(
        title = "Background Recovery",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            id = ns('sidebar'),
            width = 350,
            position = "left",
            class = "bg-light",

            # SIDEBAR -------------------------------------------------------

            bslib::accordion(
              id = "introgression_accordion",
              open = "files", # Only the first panel is open by default

              # Upload genotype file
              bslib::accordion_panel(
                title = section_title("upload", "Data Acquisition"),
                value = "files",

                fileInput(
                  inputId = ns("data_id"),
                  label = "Upload Genotype Data ",
                  placeholder = ".csv,  .xlsx",
                  multiple = FALSE,
                  accept = c(".csv", ".xlsx")
                ),
                example_data_link(ns("load_example")),
                uiOutput(ns("example_badge")),
                tags$details(
                  class = "small text-muted mb-3",
                  tags$summary("Data format"),
                  tags$ul(
                    class = "ps-3 mt-2 mb-0",
                    tags$li("One row per sample. Include both parents and the progenies."),
                    tags$li("Sample information columns first (e.g. sample ID, batch), then one column per marker, starting at Marker Start."),
                    tags$li("Marker columns hold the genotype calls, e.g. A, G or A / G for Agriplex.")
                  )
                ),

                layout_columns(
                  col_widths = c(6, 6),

                  selectInput(
                    inputId = ns("data_type"),
                    label = "Data Type",
                    choices = c("Agriplex", "Kasp / DArTag"),
                    selected = "Agriplex",
                    width = "100%"
                  ),
                  numericInput(
                    inputId = ns("marker_start"),
                    label = label_help(
                      "Marker Start",
                      "Column number where the marker (genotype call) columns begin. Columns before it are sample information."
                    ),
                    value = 7,
                    min = 1,
                    width = "100%"
                  )
                ),
                textInput(
                  inputId = ns("allele_sep"),
                  label = label_help(
                    "Allele Separator",
                    "Characters between the two alleles in a call, e.g. ' / ' in A / G. Spaces count."
                  ),
                  value = " / ",
                  width = "100%"
                )
              ),

              # Map file, or generate one
              bslib::accordion_panel(
                title = section_title("map", "Genomic Mapping Setup"),
                value = "mapfile",

                radioButtons(
                  inputId = ns("choice"),
                  label = label_help(
                    "Do you have a map file?",
                    "A map file gives each marker's chromosome and position. Without one, positions are read from the marker names (e.g. S1_12345)."
                  ),
                  choices = c("Yes" = "yes", "Generate one" = "no"),
                  selected = "yes"
                ),

                # Upload slot when they have one
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
                      label = "SNP ID Column",
                      choices = NULL,
                      width = "100%"
                    )
                  )
                ),

                conditionalPanel(
                  condition = paste0('input["', ns("choice"), '"] == "no"'),
                  layout_columns(
                    col_widths = c(6, 6),
                    textInput(
                      inputId = ns("prefix_marker"),
                      label = label_help(
                        "Marker Prefix",
                        "Text before the chromosome number in marker names, e.g. 'S' in S1_12345."
                      ),
                      value = "S",
                      width = "100%"
                    ),
                    textInput(
                      inputId = ns("sep_marker"),
                      label = label_help(
                        "Name Separator",
                        "Character between chromosome and position in marker names, e.g. '_' in S1_12345."
                      ),
                      value = "_",
                      width = "100%"
                    )
                  )
                )
              ),

              # Parents and batch
              bslib::accordion_panel(
                title = section_title("sliders", "Selection Parameters"),
                value = "settings",

                selectizeInput(
                  inputId = ns("sample_id"),
                  label = "Sample ID Column",
                  choices = NULL,
                  width = "100%"
                ),

                selectizeInput(
                  inputId = ns("cluster_by"),
                  label = label_help(
                    "Group By",
                    "Optional column used to keep only one group of samples, e.g. one batch or family."
                  ),
                  choices = NULL,
                  width = "100%"
                ),

                uiOutput(outputId = ns("cluster_focus_ui")),

                bslib::layout_columns(
                  col_widths = c(6, 6),
                  selectizeInput(
                    inputId = ns("rp"),
                    label = label_help("Recurrent Parent", "The elite parent the progenies are backcrossed to."),
                    choices = NULL,
                    width = "100%"
                  ),
                  selectizeInput(
                    inputId = ns("dp"),
                    label = label_help("Donor Parent", "The parent that carries the trait being introgressed."),
                    choices = NULL,
                    width = "100%"
                  )
                )
              ),

              # Data cleaning
              bslib::accordion_panel(
                title = section_title("filter", "Genotypic Data Cleaning"),
                value = "qc",

                bslib::input_switch(
                  id = ns("apply_par_poly"),
                  label = label_help(
                    "Remove Monomorphic Parents",
                    "Drop markers where both parents have the same genotype; they cannot tell the parents apart."
                  ),
                  value = TRUE
                ),

                bslib::input_switch(
                  id = ns("apply_par_miss"),
                  label = label_help("Remove Missing Parent Data", "Drop markers where either parent has no call."),
                  value = TRUE
                ),

                bslib::input_switch(
                  id = ns("apply_geno_good"),
                  label = label_help(
                    "Apply Genotype Error Check",
                    "Flag progeny calls that are impossible given the parents' genotypes."
                  ),
                  value = TRUE
                ),

                bslib::input_switch(
                  id = ns("apply_par_homo"),
                  label = label_help("Filter Heterozygous Parents", "Drop markers where either parent is heterozygous."),
                  value = TRUE
                )
              )
            ),

            run_button(ns("config"), "Process Data", icon_name = "gears")
          ),

          # PLACEHOLDER ---------------------------------------------------

          # Shown until the data has been processed
          div(
            id = ns("placeholder_ui"),
            class = "text-center text-muted py-5",
            div(class = "h1 mb-3", icon("chart-bar")),
            tags$p(
              "Upload backcross genotype data in the sidebar, or load the example data, then click",
              tags$b("Process Data"), "."
            ),
            tags$p(
              class = "small",
              "You'll then see how much of the recurrent parent's genome each progeny has recovered (RPP)."
            )
          ),

          # RESULTS -------------------------------------------------------

          # Hidden until the data has been processed
          shinyjs::hidden(
            div(
              id = ns("analysis_ui"),
              bslib::navset_card_underline(
                id = ns("config_pages"),

                # RPP settings tab
                bslib::nav_panel(
                  title = "RPP Settings",
                  value = "main_config",
                  icon = icon("sliders"),
                  bslib::layout_columns(
                    col_widths = bslib::breakpoints(sm = 12, lg = c(6, 6)),

                    div(
                      tags$p(class = "fw-semibold mb-2", "Map columns"),
                      selectInput(inputId = ns("snp_ids"), label = "Marker Column", choices = NULL, width = "100%"),
                      selectInput(inputId = ns("chr"), label = "Chromosome Column", choices = NULL, width = "100%"),
                      selectInput(inputId = ns("chr_pos"), label = "Position Column", choices = NULL, width = "100%"),
                      radioButtons(
                        inputId = ns("weight_rpp"),
                        label = label_help(
                          "Weighted RPP",
                          "Yes: each marker counts in proportion to the distance it covers, so unevenly spaced markers don't skew the result. No: every marker counts equally."
                        ),
                        choices = c("Yes" = TRUE, "No" = FALSE),
                        selected = TRUE,
                        inline = TRUE
                      )
                    ),

                    div(
                      tags$p(class = "fw-semibold mb-2", "Selection threshold"),
                      radioButtons(
                        inputId = ns("rpp_threshold_mode"),
                        label = label_help(
                          "Threshold based on",
                          "Fixed: a value you choose. BC generation: the RPP expected on average for that backcross generation (e.g. 93.75% for BC3)."
                        ),
                        choices = c("Fixed value" = "fixed", "BC generation" = "bc_gen"),
                        selected = "fixed",
                        inline = TRUE
                      ),
                      conditionalPanel(
                        condition = paste0('input["', ns("rpp_threshold_mode"), '"] == "fixed"'),
                        numericInput(
                          inputId = ns("rpp_threshold"),
                          label = "Minimum RPP (0-1)",
                          value = 0.93, min = 0, max = 1, step = 0.01, width = "100%"
                        )
                      ),
                      conditionalPanel(
                        condition = paste0('input["', ns("rpp_threshold_mode"), '"] == "bc_gen"'),
                        numericInput(inputId = ns("bc_gen"), label = "BC Generation", value = 3, min = 1, width = "100%")
                      ),
                      bslib::input_switch(
                        id = ns("show_above_thresh"),
                        label = "Plot only progenies above the threshold",
                        value = FALSE
                      ),
                      selectizeInput(
                        inputId = ns("rpp_replicates"),
                        label = label_help(
                          "Parent replicates",
                          "Extra samples of the parents (e.g. BTx623b). They stay in the RPP table as controls but are left out of the progeny counts."
                        ),
                        choices = NULL,
                        multiple = TRUE,
                        width = "100%"
                      ),
                      bslib::accordion(
                        open = FALSE,
                        class = "mt-2",
                        bslib::accordion_panel(
                          title = "Result columns",
                          icon = icon("table-columns"),
                          tags$p(class = "small text-muted", "Picked automatically; change only if needed."),
                          selectInput(inputId = ns("rpp_col"), label = "Total RPP Column", choices = NULL, width = "100%"),
                          selectInput(inputId = ns("rpp_sample_id"), label = "Sample ID Column", choices = NULL, width = "100%")
                        )
                      )
                    )
                  ),
                  div(
                    class = "d-flex justify-content-end mt-2",
                    div(
                      style = "min-width: 260px;",
                      run_button(ns("generate_rpp_plot"), "Generate Plot")
                    )
                  )
                ),

                # RPP results tab: plot, then table
                bslib::nav_panel(
                  title = "RPP Results",
                  value = "plot_tab",
                  icon = icon("chart-bar"),
                  ui_result_card(
                    "Recurrent parent percentage (RPP) per progeny",
                    full_screen = TRUE,
                    actions = tagList(
                      plot_options_button(
                        bslib::layout_columns(
                          colourpicker::colourInput(inputId = ns("bar_col"), label = "Bar Colour", value = "#6495ED"),
                          colourpicker::colourInput(inputId = ns("thresh_line_col"), label = "Threshold Line", value = "#B22222")
                        ),
                        bslib::layout_columns(
                          numericInput(inputId = ns("alpha"), label = "Opacity", value = 0.9, min = 0, max = 1, step = 0.1),
                          numericInput(inputId = ns("text_size"), label = "Text Size", value = 15, min = 1)
                        ),
                        bslib::layout_columns(
                          numericInput(inputId = ns("bar_width"), label = "Bar Width", value = 0.5, min = 0.1, step = 0.1),
                          numericInput(inputId = ns("aspect_ratio"), label = "Aspect Ratio", value = 0.5, min = 0.1, step = 0.1)
                        ),
                        numericInput(inputId = ns("text_scale_fct"), label = "Label Text Scaling", value = 0.1, min = 0.1, step = 0.1),
                        actionButton(
                          inputId = ns("redraw_rpp"),
                          label = "Apply",
                          icon = icon("rotate"),
                          class = "btn-success w-100"
                        ),
                        title = "Plot appearance"
                      ),
                      plot_download_button(
                        ns, "download_plot2",
                        file_name = "rpp_barplot",
                        formats = c("PDF" = "pdf"),
                        width = 8, height = 6
                      )
                    ),
                    shinycssloaders::withSpinner(
                      plotOutput(outputId = ns("rpp_bar"), width = "100%", height = "600px"),
                      type = 4, color = "#0dc5c1"
                    )
                  ),
                  ui_result_card(
                    "RPP table",
                    full_screen = TRUE,
                    actions = downloadButton(ns("download_rpp_table"), "CSV", class = "btn-sm btn-success"),
                    uiOutput(ns("rpp_summary")),
                    DT::DTOutput(ns("rpp_table"))
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

    # MAP FILE GUIDE --------------------------------------------------------

    # Set when the example turns on "Generate one", so the guide is skipped
    skip_choice_modal <- FALSE

    # Marker naming guide when they choose "Generate one"
    observe({
      if (input$choice == "no") {
        # Don't show the naming guide when the example data switched this on
        if (skip_choice_modal) {
          skip_choice_modal <<- FALSE
          return()
        }
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

              # Naming rules
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

              # Example names
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

    # DATA UPLOAD -----------------------------------------------------------

    # Genotype file to read: the latest upload, or the bundled example
    data_path <- reactiveVal(NULL)
    example_active <- reactiveVal(FALSE)
    example_name <- "agriplex_dat.csv"

    # New upload
    observeEvent(input$data_id, {
      data_path(input$data_id$datapath)
      example_active(FALSE)
    })

    # Read the genotype file
    validated_data <- reactive({
      req(data_path())
      read_mapfile(filepath = data_path())
    })

    # EXAMPLE DATA ----------------------------------------------------------

    # Vignette Agriplex example: stg5 NILs in Batch 3, BTx623 x BTx642,
    # map file generated from the marker names
    observeEvent(input$load_example, {
      updateSelectInput(session, "data_type", selected = "Agriplex")
      updateNumericInput(session, "marker_start", value = 7)
      updateTextInput(session, "allele_sep", value = " / ")
      if (!identical(input$choice, "no")) skip_choice_modal <<- TRUE
      updateRadioButtons(session, "choice", selected = "no")
      updateTextInput(session, "prefix_marker", value = "S")
      updateTextInput(session, "sep_marker", value = "_")

      data_path(example_file(example_name))
      example_active(TRUE)
    })

    # Clear the example
    observeEvent(input$clear_example, {
      data_path(NULL)
      example_active(FALSE)
    })

    # "Using example data" badge
    output$example_badge <- renderUI({
      req(example_active())
      example_data_badge(
        example_name,
        clear_id = ns("clear_example"),
        note = "Batch 3 (stg5 NILs) is selected with BTx623a and BTx642a as parents. Click Process Data to run."
      )
    })

    # COLUMN DROPDOWNS AND PARENTS -------------------------------------------

    # Column names of the genotype file
    data_colnames <- reactive({
      req(validated_data())
      colnames(validated_data())
    })

    # Fill the batch and genotype column dropdowns
    observe({
      req(data_colnames())

      # Batch column
      updateSelectizeInput(
        session,
        inputId = "cluster_by",
        choices = c("None", data_colnames()),
        selected = if (isTRUE(example_active())) "Batch" else NULL,
        server = TRUE
      )

      # Genotype (sample name) column
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

    # Fill the parent pickers with sample names, mapped to their row index
    observe({
      req(sub_validated_data(), input$sample_id)

      sample_names <- sub_validated_data()[[input$sample_id]]
      parent_choices <- stats::setNames(as.character(seq_along(sample_names)), sample_names)

      updateSelectizeInput(session,
        inputId = "rp",
        choices = parent_choices,
        selected = parent_choices[1]
      )

      updateSelectizeInput(session,
        inputId = "dp",
        choices = parent_choices,
        selected = if (length(parent_choices) >= 3) parent_choices[3] else parent_choices[1]
      )
    })

    # BATCH FILTER ----------------------------------------------------------

    # Batch value picker, shown when a batch column is chosen
    output$cluster_focus_ui <- renderUI({
      # Nothing to pick when no batch column is chosen
      req(input$cluster_by, validated_data())
      if (input$cluster_by == "None") {
        return(NULL)
      }

      # Batch values
      choices_to_focus <- unique(validated_data()[[input$cluster_by]])

      selectInput(
        inputId = ns("cluster_focus"),
        label = paste("Select", input$cluster_by, "to focus on:"),
        choices = choices_to_focus,
        selected = if (isTRUE(isolate(example_active())) && "3" %in% as.character(choices_to_focus)) "3" else NULL
      )
    })
    # Render even while its accordion panel is collapsed
    outputOptions(output, "cluster_focus_ui", suspendWhenHidden = FALSE)

    # Keep only the chosen batch
    sub_validated_data <- reactive({
      req(validated_data(), input$cluster_by)

      if (!input$cluster_by %in% colnames(validated_data()) || input$cluster_by == "None") {
        return(validated_data())
      } else {
        req(input$cluster_focus)
        return(validated_data()[validated_data()[[input$cluster_by]] == input$cluster_focus, ])
      }
    })

    # MAP FILE --------------------------------------------------------------

    # Read map file if user has one
    map_file <- reactive({
      req(input$mapfile)
      read_mapfile(filepath = input$mapfile$datapath)
    })

    # Fill the SNP ID column dropdown from the map file
    observe({
      req(map_file())
      updateSelectizeInput(session,
        inputId = "snp_id", server = TRUE,
        choices = colnames(map_file()),
        selected = safe_grep_match(
          pattern = "id",
          choices = colnames(map_file())
        )
      )
    })

    # PROCESS DATA ----------------------------------------------------------

    # Check required fields first, so a missing one shows a modal naming it
    observeEvent(input$config, {
      if (!validate_required_inputs(list(
        "Uploaded Genotype Data (upload a file first)" = sub_validated_data(),
        "Marker Start" = input$marker_start,
        "Sample ID Column" = input$sample_id,
        "Allele Separator" = input$allele_sep,
        "Recurrent Parent" = input$rp,
        "Donor Parent" = input$dp,
        "Map File Source" = input$choice
      ))) {
        return()
      }
    })

    # Process data on Process Data
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
          rp = as.numeric(input$rp),
          dp = as.numeric(input$dp),
          Prefix = if (input$choice == "no") input$prefix_marker else NULL,
          feedback = input$choice,
          mapfile_path = if (!is.null(map_file())) map_file() else NULL
        )
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Couldn't process the data",
          text = paste0(e$message, " Check the data type, marker start column and parents in the sidebar."),
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

    # Show the results and tell them what to do next
    observeEvent(Result(), {
      req(Result()) # Only proceed if Result is not NULL

      # Swap placeholder for results
      shinyjs::hide("placeholder_ui")
      shinyjs::show("analysis_ui")

      # Close the sidebar for more room
      bslib::toggle_sidebar(id = 'sidebar' , open = "closed")

      shinyWidgets::show_alert(
        title = "Data processed",
        text = "Choose the map columns and threshold under RPP Settings, then click Generate Plot.",
        type = "success",
        btn_labels = "Got it!",
        closeOnClickOutside = FALSE,
        showCloseButton = TRUE
      )
    })

    # RPP SETTINGS ----------------------------------------------------------

    # Fill the map column dropdowns
    observe({
      req(Result())

      # SNP ID
      updateSelectInput(session,
        inputId = "snp_ids",
        choices = colnames(Result()$mapfile),
        selected = safe_grep_match(
          pattern = "id",
          choices = colnames(Result()$mapfile)
        )
      )

      # Chromosome
      updateSelectInput(session,
        inputId = "chr",
        choices = colnames(Result()$mapfile),
        selected = safe_grep_match(
          pattern = "chr",
          choices = colnames(Result()$mapfile)
        )
      )

      # Position
      updateSelectInput(session,
        inputId = "chr_pos",
        choices = colnames(Result()$mapfile),
        selected = safe_grep_match(
          pattern = "pos",
          choices = colnames(Result()$mapfile)
        )
      )
    })

    # RPP CALCULATION -------------------------------------------------------

    # Run calc_rpp_bc()
    calc_rpp_bc_result <- reactive({
      req(
        Result(), input$chr_pos, input$chr, input$weight_rpp,
        input$snp_ids, sub_validated_data(), input$rp, input$sample_id
      )

      # Recurrent parent name from its row index
      parent_name <- sub_validated_data()[as.numeric(input$rp), input$sample_id]

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

    # RPP PLOT --------------------------------------------------------------

    # Fill the RPP plot dropdowns
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

    # Go to the results tab on Generate Plot
    observeEvent(input$generate_rpp_plot, {
      if (!validate_required_inputs(list(
        "RPP Results (calculate RPP first)" = calc_rpp_bc_result(),
        "Total RPP Column" = input$rpp_col,
        "Sample ID Column" = input$rpp_sample_id,
        "Selection Threshold Basis" = input$rpp_threshold_mode
      ))) {
        return()
      }
      bslib::nav_select("config_pages", "plot_tab")
    })

    # Bar plot, redrawn by Generate Plot or by Apply in Appearance.
    # No ignoreInit: the first run happens after the first click, and
    # ignoreInit would throw that click away
    rpp_barplot_result <- eventReactive(list(input$generate_rpp_plot, input$redraw_rpp), {
      req(input$generate_rpp_plot > 0)
      req(calc_rpp_bc_result(), input$rpp_col, input$rpp_sample_id, input$rpp_threshold_mode)

      tryCatch(
        {
          rpp_barplot(
            rpp_df = calc_rpp_bc_result(),
            rpp_sample_id = input$rpp_sample_id,
            rpp_col = input$rpp_col,
            rpp_threshold = if (input$rpp_threshold_mode == "fixed") input$rpp_threshold else NULL,
            text_size = input$text_size,
            text_scale_fct = input$text_scale_fct,
            alpha = input$alpha,
            bar_width = input$bar_width,
            aspect_ratio = input$aspect_ratio,
            bar_col = input$bar_col,
            thresh_line_col = input$thresh_line_col,
            show_above_thresh = input$show_above_thresh,
            bc_gen = if (input$rpp_threshold_mode == "bc_gen") input$bc_gen else NULL,
            pdf = FALSE
          )
        },
        error = function(e) {
          showNotification(paste("Couldn't draw the RPP plot:", e$message), type = "error")
          return(NULL)
        }
      )
    })

    # Show the plot
    output$rpp_bar <- renderPlot({
      tryCatch(
        {
          plot_obj <- rpp_barplot_result()
          req(plot_obj)

          if (inherits(plot_obj, "ggplot")) {
            print(plot_obj)
          } else {
            plot_obj
          }
        },
        error = function(e) {
          # req() signals a silent error; let Shiny handle it instead of alerting
          if (inherits(e, "shiny.silent.error")) stop(e)
          shinyWidgets::show_alert(
            title = "Couldn't draw the plot",
            text = e$message,
            type = "error"
          )
          return(NULL)
        }
      )
    })

    # PARENT REPLICATES -----------------------------------------------------

    # Any non-parent sample can be a replicate. Preselect names that differ
    # from a parent only by a trailing letter (e.g. BTx623a / BTx623b)
    observeEvent(Result(), {
      parents <- c(Result()$rp_name, Result()$dp_name)
      others <- setdiff(rownames(Result()$proc_kasp_f), parents)
      lettered <- parents[grepl("[A-Za-z]$", parents)]
      stem <- function(v) sub("[A-Za-z]$", "", v)
      guess <- others[grepl("[A-Za-z]$", others) & stem(others) %in% stem(lettered)]
      updateSelectizeInput(session, "rpp_replicates", choices = others, selected = guess)
    })

    # RPP TABLE -------------------------------------------------------------

    # Threshold: typed in, or the expected RPP for the BC generation
    rpp_threshold_value <- reactive({
      if (identical(input$rpp_threshold_mode, "bc_gen")) {
        req(input$bc_gen)
        unname(calc_rpp_exp(bc_gen = input$bc_gen))
      } else {
        req(input$rpp_threshold)
        input$rpp_threshold
      }
    })

    # Every progeny with its RPP and whether it clears the threshold
    rpp_table_data <- reactive({
      df <- calc_rpp_bc_result()
      req(df, input$rpp_col %in% colnames(df), input$rpp_sample_id %in% colnames(df))
      parents <- c(Result()$rp_name, Result()$dp_name)
      sid <- input$rpp_sample_id
      df <- df[!df[[sid]] %in% parents, , drop = FALSE]
      is_rep <- df[[sid]] %in% input$rpp_replicates
      df$role <- ifelse(is_rep, "parent replicate", "progeny")
      df$status <- ifelse(is_rep, "Control",
                          ifelse(df[[input$rpp_col]] >= rpp_threshold_value(), "Above", "Below"))
      # Progenies best first, parent replicates (controls) last
      df <- df[order(is_rep, -df[[input$rpp_col]]), , drop = FALSE]
      # Key columns first: sample, status, total RPP, then per-chromosome RPP
      first <- c(sid, "status", input$rpp_col)
      df[, c(first, setdiff(colnames(df), c(first, "role")), "role"), drop = FALSE]
    })

    # Summary line above the table
    output$rpp_summary <- renderUI({
      df <- rpp_table_data()
      thr <- rpp_threshold_value()
      prog <- df[df$role == "progeny", , drop = FALSE]
      tags$p(
        class = "small text-muted",
        tags$b(sum(prog$status == "Above")), " of ", nrow(prog), " progenies reach the threshold of ",
        tags$b(paste0(round(100 * thr, 2), "%")),
        if (identical(input$rpp_threshold_mode, "bc_gen")) paste0(" (expected RPP for BC", input$bc_gen, ")"),
        ". Values are the share of each chromosome, and of the whole genome, that matches the recurrent parent.",
        if (any(df$role == "parent replicate")) " Parent replicates are listed last as controls and not counted."
      )
    })

    # RPP table, coloured by status
    output$rpp_table <- DT::renderDT({
      df <- rpp_table_data()
      shown <- df[, setdiff(colnames(df), "role"), drop = FALSE]
      num_cols <- colnames(shown)[vapply(shown, is.numeric, logical(1))]
      DT::datatable(shown, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE)) |>
        DT::formatPercentage(num_cols, digits = 1) |>
        DT::formatStyle(
          "status",
          color = DT::styleEqual(c("Above", "Below", "Control"), c("#1e7e34", "#c0392b", "#6c757d")),
          fontWeight = "bold"
        )
    })

    # DOWNLOADS -------------------------------------------------------------

    # RPP table as csv
    output$download_rpp_table <- downloadHandler(
      filename = function() "rpp_table.csv",
      content = function(file) {
        utils::write.csv(rpp_table_data(), file, row.names = FALSE)
      }
    )

    # RPP plot
    output$download_plot2 <- downloadHandler(
      filename = function() {
        req(input$download_plot2_name)
        paste0(input$download_plot2_name, ".pdf")
      },
      content = function(file) {
        req(rpp_barplot_result())

        tryCatch(
          {
            grDevices::pdf(file,
              width = input$download_plot2_width,
              height = input$download_plot2_height,
              onefile = FALSE
            )

            # Single plot
            print(rpp_barplot_result())

            grDevices::dev.off()
          },
          error = function(e) {
            shinyWidgets::show_toast(
              title = "Download failed",
              type = "error",
              text = e$message,
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
