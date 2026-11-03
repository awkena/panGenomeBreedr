#' Trait Introgression Hypothesis Testing (UI)
#'
#' @description A Shiny module UI for testing hypotheses about trait introgression
#' using marker data and visualization tools.
#'
#' @param id Unique namespace ID for the module
#'
#' @return A shiny::tagList containing the UI elements
#'
#' @importFrom shiny NS tagList fileInput radioButtons selectInput textInput
#' @importFrom shiny actionButton uiOutput icon plotOutput
#' @importFrom bslib navset_card_underline nav_panel layout_sidebar sidebar
#' @importFrom bslib card card_header card_body card_footer accordion
#' @importFrom bslib accordion_panel navset_card_tab
#' @noRd
#'
#'
mod_ds_trait_introg_hypothesis_test_ui <- function(id) {
  ns <- NS(id)

  # Accordion section title with an icon
  section_title <- function(icon_name, text) {
    div(icon(icon_name, class = "me-2"), tags$span(text, class = "fw-semibold"))
  }

  tagList(
    navset_card_underline(
      nav_panel(
        title = "Trait Introgression Hypothesis Testing",
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
            div(class = "h1 mb-3", icon("table-cells")),
            tags$p(
              "Upload genotype data in the sidebar, or load the example data, then click",
              tags$b("Process Data"), "."
            ),
            tags$p(class = "small", "Heatmap settings appear here once the data is processed.")
          ),

          # RESULTS -------------------------------------------------------

          # Hidden until the data has been processed
          shinyjs::hidden(
            div(
              id = ns("analysis_ui"),
              bslib::navset_card_underline(
                id = ns('config_pages'),

                # Heatmap settings tab
                bslib::nav_panel(
                  title = "Heatmap Settings",
                  value = "main_config",
                  icon = icon("sliders"),
                  bslib::layout_columns(
                    col_widths = bslib::breakpoints(sm = 12, lg = c(6, 6)),

                    # Columns from the map file
                    div(
                      tags$p(class = "fw-semibold mb-2", "Map columns"),
                      selectInput(
                        inputId = ns("snp_ids"),
                        label = "Marker Column",
                        choices = NULL,
                        width = "100%"
                      ),
                      selectInput(
                        inputId = ns("chr"),
                        label = "Chromosome Column",
                        choices = NULL,
                        width = "100%"
                      ),
                      selectInput(
                        inputId = ns("chr_pos"),
                        label = "Position Column",
                        choices = NULL,
                        width = "100%"
                      )
                    ),

                    # What to show
                    div(
                      tags$p(class = "fw-semibold mb-2", "Samples and annotations"),
                      selectInput(
                        inputId = ns("parents"),
                        label = label_help("Parents", "Recurrent and donor parent, shown as the top rows of the heatmap."),
                        choices = NULL,
                        width = "100%",
                        multiple = TRUE
                      ),
                      numericInput(
                        inputId = ns("group_sz"),
                        label = label_help(
                          "Progeny Batch Size",
                          "Number of progenies per heatmap. Use 0 to show all progenies in one heatmap."
                        ),
                        value = 0,
                        min = 0,
                        step = 1,
                        width = "100%"
                      ),
                      radioButtons(
                        inputId = ns("options"),
                        label = label_help(
                          "Trait annotations",
                          "Mark the positions of known trait loci (QTLs) on the heatmap."
                        ),
                        choices = c(
                          "Hide" = "no",
                          "Show" = "yes"
                        ),
                        selected = "no",
                        inline = TRUE
                      ),
                      conditionalPanel(
                        condition = "input.options == 'yes'",
                        ns = ns,
                        actionButton(
                          inputId = ns("newBtn"),
                          label = "Define QTL Coordinates",
                          icon = icon("location-crosshairs"),
                          class = "btn-outline-primary mb-3 w-100"
                        ),
                        conditionalPanel(
                          condition = "output.qtls_defined === true",
                          ns = ns,
                          radioButtons(
                            inputId = ns("focus_qtl"),
                            label = "Which QTLs to annotate?",
                            choices = c(
                              "Specific QTLs (zoom in)" = "yes",
                              "All QTLs" = "no"
                            ),
                            selected = "no"
                          ),
                          conditionalPanel(
                            condition = "input.focus_qtl == 'yes'",
                            ns = ns,
                            selectInput(
                              inputId = ns("qtls"),
                              label = "Select defined QTL(s)",
                              choices = NULL,
                              multiple = TRUE
                            )
                          )
                        )
                      )
                    )
                  ),
                  div(
                    class = "d-flex justify-content-end mt-2",
                    div(
                      style = "min-width: 260px;",
                      run_button(ns("generate_heatmap"), "Generate Heatmap")
                    )
                  )
                ),

                # QTL coordinates form (hidden tab, opened from settings or the check)
                bslib::nav_panel_hidden(
                  value = "coord_entry",
                  bslib::card(
                    class = "border-0 shadow-sm",
                    full_screen = TRUE,
                    bslib::card_header(
                      div(
                        class = "d-flex justify-content-between align-items-center",
                        tags$b("Genomic Coordinate Definition"),
                        actionButton(
                          ns("addTrait"),
                          "Add New Trait",
                          icon = icon("plus"),
                          class = "btn-success btn-sm"
                        )
                      )
                    ),
                    bslib::card_body(
                      p(
                        class = "text-muted",
                        "Define specific points or genomic ranges for heatmap annotation."
                      ),
                      div(
                        id = ns("inputs_container"),
                        makeTraitInput(1, ns) # Initial row
                      )
                    ),
                    bslib::card_footer(
                      class = "bg-white",
                      div(
                        class = "d-flex justify-content-end gap-2",
                        actionButton(
                          ns("cancel_coords"),
                          "Cancel",
                          class = "btn-light"
                        ),
                        actionButton(
                          ns("submit_coords"),
                          "Save & Apply",
                          icon = icon("check"),
                          class = "btn-primary"
                        )
                      )
                    )
                  )
                ),

                # Heatmap viewer tab
                bslib::nav_panel(
                  title = "Heatmap Viewer",
                  value = "plot_tab",
                  icon = icon("chart-area"),
                  ui_result_card(
                    "Genotype heatmap",
                    full_screen = TRUE,
                    actions = tagList(
                      div(
                        style = "min-width: 220px; margin-bottom: -1rem;",
                        selectInput(
                          inputId = ns("batch_no"),
                          label = NULL,
                          choices = NULL,
                          width = "100%"
                        )
                      ),
                      plot_options_button(
                        textInput(
                          inputId = ns("legend_title"),
                          label = "Legend Title",
                          value = "Heatmap Key",
                          width = "100%"
                        ),
                        selectInput(
                          inputId = ns("col_mapping"),
                          label = "Heatmap Colours (in legend order; empty = default)",
                          choices = grDevices::colors(),
                          multiple = TRUE,
                          width = "100%"
                        ),
                        selectInput(
                          inputId = ns("col_labels"),
                          label = "Genotype Labels (in legend order; empty = default)",
                          choices = NULL,
                          multiple = TRUE,
                          width = "100%"
                        ),
                        bslib::layout_columns(
                          colourpicker::colourInput(
                            inputId = ns("panel_fill"),
                            label = "Panel Fill",
                            value = "#CCCCCC"
                          ),
                          colourpicker::colourInput(
                            inputId = ns("panel_col"),
                            label = "Panel Border",
                            value = "#FFFFFF"
                          )
                        ),
                        bslib::layout_columns(
                          numericInput(
                            inputId = ns("alpha"),
                            label = "Opacity",
                            value = 1, min = 0, max = 1, step = 0.1
                          ),
                          numericInput(
                            inputId = ns("text_size"),
                            label = "Text Size",
                            value = 12, min = 1, step = 0.1
                          )
                        ),
                        bslib::layout_columns(
                          numericInput(
                            inputId = ns("text_scale_fct"),
                            label = "Text Scaling",
                            value = 0.3, min = 0.1, step = 0.01, max = 1
                          ),
                          numericInput(
                            inputId = ns("label_offset"),
                            label = "Trait Label Offset",
                            value = 0.4, min = -10, max = 10, step = 0.1
                          )
                        ),
                        actionButton(
                          inputId = ns("redraw_heatmap"),
                          label = "Apply",
                          icon = icon("rotate"),
                          class = "btn-success w-100"
                        ),
                        title = "Heatmap appearance"
                      ),
                      plot_download_button(
                        ns, "download_plot1",
                        file_name = "Heatmap_result",
                        formats = c("PDF (all batches)" = "pdf"),
                        width = 10, height = 7
                      )
                    ),
                    shinycssloaders::withSpinner(
                      plotOutput(
                        outputId = ns("ant_heatmap"),
                        width = "100%",
                        height = "650px"
                      ),
                      type = 4,
                      color = "#0dc5c1"
                    ),
                    tags$p(
                      class = "small text-muted mt-2 mb-0",
                      icon("circle-info"),
                      "Each cell is one marker in one sample, coloured by which parent the call matches.",
                      HTML(paste(
                        "<b>Heterozygous:</b> one allele from each parent.",
                        "<b>Monomorphic:</b> the parents have the same call, so the marker is not informative.",
                        "<b>Error:</b> a call neither parent could produce.",
                        "<b>Missing:</b> no call."
                      ))
                    )
                  )
                ),

                # Introgression check tab
                bslib::nav_panel(
                  title = "Introgression Check",
                  value = "check_tab",
                  icon = icon("clipboard-check"),
                  div(
                    class = "d-flex flex-wrap align-items-end gap-3 mb-2",
                    div(
                      style = "min-width: 240px;",
                      selectInput(
                        inputId = ns("check_expected"),
                        label = label_help(
                          "Expected genotype at the QTL",
                          "What a successful progeny should carry at the target locus."
                        ),
                        choices = c(
                          "Donor homozygous (e.g. NILs)" = "donor",
                          "Heterozygous (e.g. BC F1)" = "het",
                          "Donor or heterozygous" = "carrier"
                        ),
                        width = "100%"
                      )
                    ),
                    div(
                      style = "width: 150px;",
                      numericInput(
                        inputId = ns("check_threshold"),
                        label = label_help(
                          "Agreement (%)",
                          "Share of informative markers in the QTL that must agree before a progeny is called."
                        ),
                        value = 80, min = 1, max = 100, step = 5, width = "100%"
                      )
                    ),
                    div(
                      style = "width: 150px;",
                      numericInput(
                        inputId = ns("check_nearest"),
                        label = label_help(
                          "Nearest markers",
                          "Markers used for a point QTL, or for a QTL range that contains no markers."
                        ),
                        value = 3, min = 1, max = 50, step = 1, width = "100%"
                      )
                    ),
                    div(
                      style = "min-width: 260px; flex: 1;",
                      selectizeInput(
                        inputId = ns("check_replicates"),
                        label = label_help(
                          "Parent replicates",
                          "Extra samples of the parents (e.g. BTx623b). They stay in the table as controls but are left out of the pass/fail counts."
                        ),
                        choices = NULL,
                        multiple = TRUE,
                        width = "100%"
                      )
                    ),
                    div(
                      class = "mb-3",
                      downloadButton(ns("download_check"), "CSV", class = "btn-sm btn-success")
                    )
                  ),
                  uiOutput(ns("check_body"))
                )
              )
            )
          )
        )
      )
    )
  )
}

#' Trait Introgression Hypothesis Testing (Server)
#'
#' @description Server logic for trait introgression hypothesis testing module that:
#' 1. Processes marker data files
#' 2. Generates heatmaps of genotype data
#' 3. Allows annotation of trait positions
#' 4. Performs hypothesis testing
#'
#' @param id Unique namespace ID for the module
#'
#' @importFrom shiny moduleServer observeEvent req reactive renderUI observe
#' @importFrom shiny showModal modalDialog modalButton removeModal insertUI
#' @importFrom shiny removeUI updateSelectInput
#' @importFrom utils read.csv
#' @noRd
mod_ds_trait_introg_hypothesis_test_server <- function(id) {
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
                class = "text-muted mb-4"),

              # Naming rules
              tags$ul(
                class = "list-group list-group-flush mb-4",
                tags$li(class = "list-group-item border-0 ps-0",
                        icon("chevron-right", class = "text-primary me-2 small"),
                        "A common prefix (e.g., ", tags$code("S"), ")"),
                tags$li(class = "list-group-item border-0 ps-0",
                        icon("chevron-right", class = "text-primary me-2 small"),
                        "Chromosome number (e.g., ", tags$code("1"), ")"),
                tags$li(class = "list-group-item border-0 ps-0",
                        icon("chevron-right", class = "text-primary me-2 small"),
                        "A separator character (e.g., ", tags$code("_"), ")"),
                tags$li(class = "list-group-item border-0 ps-0",
                        icon("chevron-right", class = "text-primary me-2 small"),
                        "Position number (e.g., ", tags$code("101"), ")")
              ),

              # Example names
              div(
                class = "p-3 rounded-3 bg-light border",
                style = "border-left: 5px solid #0d6efd !important;",
                tags$b("Example: ", class = "text-primary"),
                tags$code("S1_101"),
                p(class = "small text-muted mt-2 mb-0",
                  "This indicates Prefix (S), Chromosome (1), and Position (101).")
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
    observeEvent(input$close_modal, { removeModal() })

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
        choices = c('None',data_colnames()),
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
        ),server = TRUE
      )
    })

    # Fill the parent pickers with sample names, mapped to their row index
    # in sub_validated_data(), the same data rp/dp index into later
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
      req(input$cluster_by , validated_data())
      if (input$cluster_by == "None") return(NULL)

      # Batch values
      choices_to_focus <- unique(validated_data()[[input$cluster_by]])

      selectInput(inputId = ns("cluster_focus"),
                  label = paste("Select", input$cluster_by, "to focus on:"),
                  choices = choices_to_focus,
                  selected = if (isTRUE(isolate(example_active())) && "3" %in% as.character(choices_to_focus)) "3" else NULL)
    })
    # Render even while its accordion panel is collapsed, so a selected group
    # (e.g. from the example data) is available when Process Data is clicked
    outputOptions(output, "cluster_focus_ui", suspendWhenHidden = FALSE)

    # Keep only the chosen batch
    sub_validated_data <- reactive({
      req(validated_data(),input$cluster_by)

      # Ignore a batch column that is not in the data
      if (!input$cluster_by %in% colnames(validated_data())){
        return(validated_data())
      } else if(input$cluster_by != 'None'){
        req(input$cluster_focus)
        validated_data()[validated_data()[input$cluster_by] == input$cluster_focus ,]
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
                        selected = safe_grep_match(pattern = 'id',
                                                    choices = colnames(map_file()))
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
      req(sub_validated_data(),input$marker_start,input$sample_id ,
          input$allele_sep ,input$rp, input$dp, input$choice)

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Processing... Please wait."
      )

      result <- tryCatch({
        proc_nd_map_func(data = sub_validated_data() ,
                          marker_start = input$marker_start ,
                          sample_id = input$sample_id,
                          marker_sep = if (input$choice == "no") input$sep_marker else NULL,
                          apply_par_poly = input$apply_par_poly,
                          apply_par_miss = input$apply_par_miss,
                          apply_geno_good = input$apply_geno_good,
                          apply_par_homo = input$apply_par_homo,
                          snp_id = if (input$choice == "yes") input$snp_id else NULL,
                          calls_sep = check_sep(input$allele_sep),
                          data_type = if(input$data_type == "Agriplex") 'agriplex' else if (input$data_type == "Kasp / DArTag") 'kasp',
                          rp = as.numeric(input$rp),
                          dp = as.numeric(input$dp),
                          Prefix = if (input$choice == "no") input$prefix_marker else NULL,
                          feedback = input$choice,
                          mapfile_path = if (!is.null(map_file())) map_file() else NULL)

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
        text = "Choose the map columns and options under Heatmap Settings, then click Generate Heatmap.",
        type = "success",
        btn_labels = "Got it!",
        closeOnClickOutside = FALSE,
        showCloseButton = TRUE
      )
    })

    # Order of the heatmap legend labels
    color_code_checker_res <- reactive({
      req(Result() , input$dp,input$rp )
      color_code_checker(Result()$proc_kasp_f ,
                          par_1 = rownames(Result()$proc_kasp_f)[as.numeric(input$rp)],
                          par_2 = rownames(Result()$proc_kasp_f)[as.numeric(input$dp)])
    })

    # HEATMAP SETTINGS ------------------------------------------------------

    # Fill the map column and parent dropdowns
    observe({
      req(Result())

      # SNP ID
      updateSelectInput(session,
                        inputId = "snp_ids",
                        choices = colnames(Result()$mapfile),
                        selected = safe_grep_match(pattern = 'id',
                                                    choices = colnames(Result()$mapfile))

      )

      # Chromosome
      updateSelectInput(session,
                        inputId = "chr",
                        choices = colnames(Result()$mapfile),
                        selected = safe_grep_match(pattern = 'chr',
                                                    choices = colnames(Result()$mapfile))
      )

      # Position
      updateSelectInput(session,
                        inputId = "chr_pos",
                        choices = colnames(Result()$mapfile),
                        selected = safe_grep_match(pattern = 'pos',
                                                    choices = colnames(Result()$mapfile))
      )

      req(color_code_checker_res())
      # Position
      updateSelectInput(session,
                        inputId = "col_labels",
                        choices = unname(color_code_checker_res()),
                        selected = NULL
      )

      # Update parents selection
      updateSelectInput(session,
                        inputId = "parents",
                        choices = rownames(Result()$proc_kasp_f),
                        selected = c(Result()$rp, Result()$dp)
      )

    })

    # QTL COORDINATES -------------------------------------------------------

    # QTLs entered in the form
    values <- reactiveValues(
      count = 1,
      removed_traits = c(),
      trait_data = NULL
    )

    # Open the form, and remember which tab to go back to
    coord_return <- reactiveVal("main_config")
    observeEvent(input$newBtn, {
      coord_return("main_config")
      bslib::nav_select("config_pages", "coord_entry")
    })

    # Add a QTL row
    observeEvent(input$addTrait, {
      values$count <- values$count + 1
      insertUI(
        selector = paste0("#", ns("inputs_container")),
        where = "beforeEnd",
        ui = makeTraitInput(values$count, ns)
      )
    })

    # Remove a QTL row
    observe({
      for (i in 2:values$count) {
        local({
          idx <- i
          observeEvent(input[[paste0("remove_", idx)]], {
            removeUI(selector = paste0("#", ns(paste0("trait_", idx))))
            values$removed_traits <- c(values$removed_traits, idx)
          }, ignoreInit = TRUE, once = TRUE)
        })
      }
    })

    # Save the QTLs
    observeEvent(input$submit_coords, {
      trait_list <- list()

      for (i in 1:values$count) {
        if (i %in% values$removed_traits) next

        name  <- trimws(input[[paste0("name_", i)]])
        chr   <- input[[paste0("chr_", i)]]
        types <- input[[paste0("type_", i)]] # Checkbox group values

        # Validation: Name must exist and at least one type selected
        if (nchar(name) == 0 || is.null(types)) next

        # Build vector based on checkboxes
        trait_vals <- c(chr = as.numeric(chr))

        if ("pos" %in% types) {
          trait_vals["pos"] <- as.numeric(input[[paste0("pos_", i)]])
        }
        if ("range" %in% types) {
          trait_vals["start"] <- as.numeric(input[[paste0("start_", i)]])
          trait_vals["end"]   <- as.numeric(input[[paste0("end_", i)]])
        }

        trait_list[[name]] <- trait_vals
      }

      if (length(trait_list) == 0) {
        shinyWidgets::show_alert("No traits to save", "Add at least one trait with a name and a position or region before saving.", "warning")
        return()
      }

      values$trait_data <- trait_list
      bslib::nav_select("config_pages", coord_return())

      shinyWidgets::show_toast("", "Trait positions saved.", type = "success")
    })

    # Cancel
    observeEvent(input$cancel_coords, {
      bslib::nav_select("config_pages", coord_return())

      values$trait_data <- NULL # nullify all existing trait data

      # Toast to confirm nothing was saved
      shinyWidgets::show_toast(
        title = "",
        text = "No changes made to the trait positions.",
        type = "info"
      )
    })

    # Whether any QTLs are saved (shows the QTL options)
    output$qtls_defined <- reactive({
      !is.null(values$trait_data) && length(values$trait_data) > 0
    })
    outputOptions(output, "qtls_defined", suspendWhenHidden = FALSE)

    # Fill the QTL picker when focusing on QTLs
    observe({
      req(values$trait_data)
      if(input$focus_qtl == 'yes'){
        updateSelectInput(session, inputId = 'qtls',choices = names(values$trait_data))
      } else {
        updateSelectInput(session, inputId = 'qtls',choices = NULL,selected = NULL)
      }
    })

    # HEATMAP ---------------------------------------------------------------

    # Go to the viewer on Generate Heatmap
    observeEvent(input$generate_heatmap, {
      if (!validate_required_inputs(list(
        "Processed Data (process your data first)" = Result(),
        "Display Annotations Option" = input$options,
        "Marker Column" = input$snp_ids,
        "Chromosome Column" = input$chr,
        "Position Column" = input$chr_pos,
        "Parents" = input$parents
      ))) {
        return()
      }
      bslib::nav_select("config_pages", "plot_tab")
    })

    # Heatmap, redrawn by Generate Heatmap or by Apply in Appearance.
    # No ignoreInit: the first run happens after the first click, and
    # ignoreInit would throw that click away
    heatmap_plot <- eventReactive(list(input$generate_heatmap, input$redraw_heatmap), {
      req(input$generate_heatmap > 0)
      req(input$options , Result(),input$snp_ids, input$chr,input$chr_pos ,
          input$parents  )

      # Default color mapping
      col_mapping <- c('-5' = 'deeppink', '-2' = 'cornflowerblue', '-1' = 'beige',
                        '0' = 'purple2', '0.5' = 'gold', '1' = 'coral2')

      col_labels <- c('-5' = "Missing", '-2' = 'Error', '-1' = "Monomorphic",
                      '0' = input$parents[2], '0.5' = "Heterozygous", '1' = input$parents[1])

      if(input$options == "no"){
        cross_qc_heatmap(
          col_mapping = if(is.null(input$col_mapping)) col_mapping else input$col_mapping,
          col_labels  = if(is.null(input$col_labels)) col_labels else input$col_labels,
          x = Result()$proc_kasp_f,
          map_file = Result()$mapfile,
          snp_ids = input$snp_ids,
          chr = input$chr,
          chr_pos = input$chr_pos,
          parents = input$parents,
          group_sz = if(input$group_sz == 0) nrow(Result()$proc_kasp_f) - 2 else input$group_sz,
          trait_pos = NULL,
          text_scale_fct = if(input$text_scale_fct <= 0) 0.1 else input$text_scale_fct,
          legend_title = if(is.null(input$legend_title)) "Genetic Heatmap Key" else input$legend_title,
          text_size = if(is.null(input$text_size)) 1 else input$text_size ,
          alpha = if(is.null(input$alpha)) 1 else input$alpha,
          panel_fill = input$panel_fill,
          panel_col = input$panel_col,
          # Suppose user forgets to set it
          label_offset = if(is.null(input$label_offset)) -1 else input$label_offset,
          pdf = FALSE
        )
      } else if(input$options == 'yes' && input$focus_qtl == 'yes' ){ # Zoomed in plot
        req(values$trait_data ,input$qtls)

        df <- Result()$proc_kasp_f
        mapfile <- Result()$mapfile

        # Automatically setting the start and end from user defined position
        trait_lists <- values$trait_data
        target_qtl <- input$qtls # user's target or focused qtls

        # Extract genomic regions of interest. based on user selected list.
        focus_trait_list <- lapply(trait_lists[target_qtl], function(x) {
          if (all(c("start", "end") %in% names(x))) {
            return(c(chr = unname(x["chr"]), start = unname(x["start"]), end = unname(x["end"])))
          }
          else if ("pos" %in% names(x)) {
            return(c(chr = unname(x["chr"]), start = unname(x["pos"]), end = unname(x["pos"])))
          }
        })

        # Get marker range
        marker_indices <- unique(unlist(lapply(focus_trait_list, function(qtl) {
          which(
            mapfile[[input$chr]] == qtl["chr"] &
              mapfile[[input$chr_pos]] >= qtl["start"] &
              mapfile[[input$chr_pos]] <= qtl["end"]
          )
        })))

        marker_indices <- unique(unlist(lapply(focus_trait_list, function(qtl) {

          # Get all markers on the correct chromosome
          chr_markers <- which(mapfile[[input$chr]] == qtl["chr"])

          # Check if there are any markers at all on this chromosome
          if (length(chr_markers) == 0) return(NULL)

          # Find markers strictly within the range
          in_range <- which(
            mapfile[[input$chr]] == qtl["chr"] &
              mapfile[[input$chr_pos]] >= qtl["start"] &
              mapfile[[input$chr_pos]] <= qtl["end"]
          )

          # if range is empty grab the nearest available
          if (length(in_range) < 3) {
            midpoint <- (qtl["start"] + qtl["end"]) / 2
            distances <- abs(mapfile[chr_markers, "pos"] - midpoint)

            # Determine how many to pull
            num_to_pull <- min(20, length(chr_markers)) #pull nearest 20 markers, if not suffcient take whats left

            # Get the indices of the closest markers without exceeding availability
            in_range <- chr_markers[order(distances)[1:num_to_pull]]
          }

          return(in_range)
        })))

        # Sort the range in order.
        marker_indices <- sort(marker_indices)

        # Subset the genotype data and the mapfile
        df_updated <- df[, marker_indices, drop = F ]
        subset_map <- mapfile[marker_indices, , drop = F]

        cross_qc_annotate(x = df_updated ,
                          map_file = subset_map,
                          snp_ids = input$snp_ids,
                          chr = input$chr,
                          chr_pos = input$chr_pos,
                          parents = input$parents,
                          trait_pos =  focus_trait_list,
                          group_sz = if(input$group_sz == 0) nrow(df_updated) - 2 else input$group_sz,
                          pdf = F,
                          legend_title = if(is.null(input$legend_title)) "Genetic Heatmap Key" else input$legend_title,
                          col_mapping = if(is.null(input$col_mapping)) col_mapping else input$col_mapping,
                          col_labels = if(is.null(input$col_labels)) col_labels else input$col_labels,
                          panel_fill = input$panel_fill ,
                          panel_col = input$panel_col,
                          alpha = input$alpha,
                          text_size = if(is.null(input$text_size)) 12 else input$text_size,
                          text_scale_fct = if(input$text_scale_fct <= 0) 0.1 else input$text_scale_fct,
                          label_offset = if(is.null(input$label_offset)) -1 else input$label_offset
        )
        }
      else if(input$options == "yes" && input$focus_qtl == 'no'){
        req(values$trait_data)

        # Extract genomic positions.
        positions_only_list <- lapply(values$trait_data, function(x) {
          if ("pos" %in% names(x)) {
            return(c(chr = unname(x["chr"]), pos = unname(x["pos"])))
          }
          else if (all(c("start", "end") %in% names(x))) {
            midpoint <- (x["start"] + x["end"]) / 2
            return(c(chr = unname(x["chr"]), pos = unname(midpoint)))
          }
        })

        # Cross qc heatmap2: a single vertical bar
        cross_qc_heatmap2(
          col_mapping = if(is.null(input$col_mapping)) col_mapping else input$col_mapping,
          col_labels  = if(is.null(input$col_labels)) col_labels else input$col_labels,
          x = Result()$proc_kasp_f,
          map_file = Result()$mapfile,
          snp_ids = input$snp_ids,
          chr = input$chr,
          chr_pos = input$chr_pos,
          parents = input$parents,
          trait_pos = positions_only_list,
          group_sz = if(input$group_sz == 0) nrow(Result()$proc_kasp_f) - 2 else input$group_sz,
          text_scale_fct = if(input$text_scale_fct <= 0) 0.1 else input$text_scale_fct,
          legend_title = if(is.null(input$legend_title)) "Genetic Heatmap Key" else input$legend_title,
          text_size = if(is.null(input$text_size)) 12 else input$text_size,
          alpha = if(is.null(input$alpha)) 1 else input$alpha,
          label_offset = if(is.null(input$label_offset)) -1 else input$label_offset,
          panel_fill = input$panel_fill,
          panel_col = input$panel_col,
          pdf = FALSE
        )
      }

    })

    # Fill the progeny batch picker
    observe({
      req(heatmap_plot())
      updateSelectInput(session,inputId = 'batch_no',choices = names(heatmap_plot()))
    })

    # Heatmap for the selected batch
    heatmap_1 <- reactive({
      req(heatmap_plot(), input$batch_no)
      if (input$batch_no %in% names(heatmap_plot())) {
        heatmap_plot()[[input$batch_no]]
      } else {
        NULL
      }
    })

    # Show the heatmap
    output$ant_heatmap <- renderPlot({
      tryCatch({
        plot_obj <- heatmap_1()
        req(plot_obj)

        if (inherits(plot_obj, "ggplot")) {
          print(plot_obj)

        } else {

          plot_obj
        }

      }, error = function(e) {
        # req() signals a silent error; let Shiny handle it instead of alerting
        if (inherits(e, "shiny.silent.error")) stop(e)
        shinyWidgets::show_alert(
          title = "Couldn't draw the heatmap",
          text = e$message,
          type = "error"
        )
        NULL
      })
    })

    # INTROGRESSION CHECK ---------------------------------------------------

    # Does each progeny carry the donor segment at the QTLs?

    # Any non-parent sample can be a parent replicate. Preselect names that
    # differ from a parent only by a trailing letter (e.g. BTx623a / BTx623b)
    observeEvent(Result(), {
      parents <- c(Result()$rp_name, Result()$dp_name)
      others <- setdiff(rownames(Result()$proc_kasp_f), parents)
      lettered <- parents[grepl("[A-Za-z]$", parents)]
      stem <- function(v) sub("[A-Za-z]$", "", v)
      guess <- others[grepl("[A-Za-z]$", others) & stem(others) %in% stem(lettered)]
      updateSelectizeInput(session, "check_replicates", choices = others, selected = guess)
    })

    # Run the check for every progeny and QTL
    check_result <- reactive({
      req(Result(), values$trait_data, input$check_expected,
          input$check_threshold, input$check_nearest)

      map <- Result()$mapfile
      cols <- colnames(map)
      pick <- function(sel, pattern) {
        if (!is.null(sel) && sel %in% cols) sel else safe_grep_match(pattern, choices = cols)[1]
      }

      tryCatch(
        introgression_check(
          x = Result()$proc_kasp_f,
          map_file = map,
          snp_ids = pick(input$snp_ids, "id"),
          chr = pick(input$chr, "chr"),
          chr_pos = pick(input$chr_pos, "pos"),
          parents = c(Result()$rp_name, Result()$dp_name),
          trait_pos = values$trait_data,
          expected = input$check_expected,
          threshold = min(max(input$check_threshold, 1), 100) / 100,
          n_nearest = max(1, round(input$check_nearest))
        ),
        error = function(e) {
          structure(list(message = conditionMessage(e)), class = "check_error")
        }
      )
    })

    # Value boxes and table, or a prompt to define QTLs first
    output$check_body <- renderUI({
      if (is.null(values$trait_data)) {
        return(div(
          class = "text-center text-muted py-5",
          div(class = "h1 mb-3", icon("location-crosshairs")),
          tags$p("Define the QTL position(s) to test whether each progeny carries the donor segment there."),
          actionButton(
            ns("check_define_qtl"), "Define QTL Coordinates",
            icon = icon("location-crosshairs"),
            class = "btn-outline-primary"
          )
        ))
      }

      res <- check_result()
      if (inherits(res, "check_error")) {
        return(div(class = "alert alert-danger", "Could not run the check: ", res$message))
      }

      res <- res[!res$sample %in% input$check_replicates, , drop = FALSE]
      n_prog <- length(unique(res$sample))
      by_sample <- tapply(res$pass, res$sample, function(p) {
        if (any(is.na(p))) NA else all(p)
      })
      n_pass <- sum(by_sample %in% TRUE)
      n_fail <- sum(by_sample %in% FALSE)
      n_na <- sum(is.na(by_sample))

      tagList(
        bslib::layout_column_wrap(
          width = "170px",
          fill = FALSE,
          bslib::value_box("Progenies tested", n_prog, showcase = icon("seedling"),
                            showcase_layout = bslib::showcase_top_right()),
          bslib::value_box("Pass at all QTLs", n_pass, showcase = icon("circle-check"),
                            showcase_layout = bslib::showcase_top_right()),
          bslib::value_box("Fail", n_fail, showcase = icon("circle-xmark"),
                            showcase_layout = bslib::showcase_top_right()),
          bslib::value_box("No informative markers", n_na, showcase = icon("circle-question"),
                            showcase_layout = bslib::showcase_top_right())
        ),
        tags$p(
          class = "small text-muted mt-2",
          "A progeny passes when its call at every QTL matches the expected genotype.",
          "Donor, Het and Recurrent are the shares of informative markers in the QTL.",
          if (length(input$check_replicates)) "Parent replicates are listed as controls and not counted."
        ),
        DT::DTOutput(ns("check_table"))
      )
    })

    # Check table
    output$check_table <- DT::renderDT({
      res <- check_result()
      req(is.data.frame(res))
      shown <- res
      shown$pass <- ifelse(is.na(shown$pass), "-", ifelse(shown$pass, "Pass", "Fail"))
      shown$pass[shown$sample %in% input$check_replicates] <- "Control"
      # Controls last
      shown <- shown[order(shown$sample %in% input$check_replicates), , drop = FALSE]
      names(shown) <- c("Sample", "QTL", "Chr", "Region (bp)", "Markers from",
                        "Markers", "Informative", "Donor", "Het", "Recurrent",
                        "Call", "Result")
      DT::datatable(shown, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE)) |>
        DT::formatPercentage(c("Donor", "Het", "Recurrent"), digits = 0) |>
        DT::formatStyle(
          "Result",
          color = DT::styleEqual(c("Pass", "Fail", "Control"), c("#1e7e34", "#c0392b", "#6c757d")),
          fontWeight = "bold"
        )
    })

    # Define QTL Coordinates button on this tab
    observeEvent(input$check_define_qtl, {
      coord_return("check_tab")
      bslib::nav_select("config_pages", "coord_entry")
    })

    # DOWNLOADS -------------------------------------------------------------

    # Check table as csv
    output$download_check <- downloadHandler(
      filename = function() "introgression_check.csv",
      content = function(file) {
        res <- check_result()
        req(is.data.frame(res))
        res$role <- ifelse(res$sample %in% input$check_replicates, "parent replicate", "progeny")
        utils::write.csv(res, file, row.names = FALSE)
      }
    )

    # Heatmaps, one page per batch
    output$download_plot1 <- downloadHandler(
      filename = function() {
        req(input$download_plot1_name)
        paste0(input$download_plot1_name, ".pdf")  # Single PDF containing ALL plots
      },
      content = function(file) {
        req(heatmap_plot())  # Ensure plots exist

        tryCatch(
          {
            # Start PDF (onefile=TRUE ensures multi-page)
            grDevices::pdf(file,
                            width = input$download_plot1_width,
                            height = input$download_plot1_height,
                            onefile = TRUE)

            # Print ALL plots regardless of type
            for (plot_obj in heatmap_plot()) {
              if (inherits(plot_obj, "grob")) {
                grid::grid.draw(plot_obj)  # Handle grid graphics
              } else {
                print(plot_obj)  # Handle ggplot2/base R plots
              }
            }

            grDevices::dev.off()
          },
          error = function(e) {
            shinyWidgets::show_toast(
              title = '' ,
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
