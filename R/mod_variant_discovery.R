#' Variant Discovery UI Function
#'
#' @description A shiny Module for variant discovery that allows database connections
#'   and various query operations for genomic variants.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @return A UI for Variant Discovery
#'
#' @noRd
#'
#' @importFrom shiny NS tagList textInput actionButton
#' @importFrom shiny uiOutput selectInput numericInput verbatimTextOutput
#' @importFrom shiny icon div h4 span helpText tags fluidRow column br
#' @importFrom bslib navset_card_underline sidebar card card_header
#' @importFrom bslib card_body card_footer nav_panel navset_card_tab
#' @importFrom bslib layout_column_wrap input_switch layout_columns navset_card_pill nav_panel_hidden
#'
mod_variant_discovery_ui <- function(id) {
  ns <- NS(id)

  #------------------------------------------------------#
  #                   Database Info Tab
  #------------------------------------------------------#
  info_tab <- function(ns) {
    tagList(
      
      # Card for Leaflet Map
      bslib::card(
        class = "shadow-sm mb-3",
        bslib::card_header(
          icon("map-marked-alt", class = "me-1 text-primary"),
          "Geographic Distribution of 1676 Sorghum Accessions"
        ),
        bslib::card_body(
          leaflet::leafletOutput(ns("accession_map"), height = "500px")
        ),
        bslib::card_footer(
          selectInput(
            ns("color_by_col"),
            "Color Points By:",
            choices = NULL, 
            width = "300px"
          )
        )
      ),
      # Card 1: Variant Impact Summary 
      bslib::card(
        class = "shadow-sm mb-3", 
        bslib::card_header(
          class = "d-flex justify-content-between align-items-center",
          tagList(
            icon("chart-pie", class = "me-1 text-primary"),
            "Variant Impact Summary"
          ),
          shinyWidgets::actionBttn(
            ns("get_impact"),
            "Get Info",
            style = "bordered",
            size = "xs",
            color = "primary",
            icon = icon("rotate")
          )
        ),
        bslib::card_body(
          reactable::reactableOutput(ns("table_impact_id"))
        )
      ),
      # Card 2: Variant Statistics 
      bslib::card(
        class = "shadow-sm mb-3",
        bslib::card_header(
          class = "d-flex justify-content-between align-items-center",
          tagList(
            icon("chart-bar", class = "me-1 text-primary"),
            "Variant Statistics"
          ),
          shinyWidgets::actionBttn(
            ns("get_stats"),
            "Get Info",
            style = "bordered",
            size = "xs",
            color = "primary",
            icon = icon("rotate")
          )
        ),
        bslib::card_body(
          reactable::reactableOutput(ns("table_var_stats_id"))
        )
      ),
      # Row with Cards 3, 4, and 5
      bslib::layout_columns(
        col_widths = c(6, 6, -2, 8, -2),
        # Summarised Tables
        bslib::card(
          class = "shadow-sm",
          bslib::card_header(
            class = "d-flex justify-content-between align-items-center",
            tagList(
              icon("table", class = "me-1 text-primary"),
              "Summarised Tables"
            ),
            shinyWidgets::actionBttn(
              ns("get_summary"),
              "Get Info",
              style = "bordered",
              size = "xs",
              color = "primary",
              icon = icon("rotate")
            )
          ),
          bslib::card_body(
            reactable::reactableOutput(ns("sum_sqlite_id"))
          )
        ),
        # Variant Type Count
        bslib::card(
          class = "shadow-sm",
          bslib::card_header(
            class = "d-flex justify-content-between align-items-center",
            tagList(
              icon("tags", class = "me-1 text-primary"),
              "Variant Type Count"
            ),
            shinyWidgets::actionBttn(
              ns("get_types"),
              "Get Info",
              style = "bordered",
              size = "xs",
              color = "primary",
              icon = icon("rotate")
            )
          ),
          bslib::card_body(
            reactable::reactableOutput(ns("count_variant_typ_id"))
          )
        ),
        # Inspect Table Schema
        bslib::card(
          class = "shadow-sm",
          bslib::card_header(
            class = "d-flex justify-content-between align-items-center",
            tagList(
              icon("list-ul", class = "me-1 text-primary"),
              "Available Data Fields"
            ),
            shinyWidgets::actionBttn(
              ns("get_schema"),
              "Get Info",
              style = "bordered",
              size = "xs",
              color = "primary",
              icon = icon("rotate")
            )
          ),
          bslib::card_body(
            selectInput(
              ns("table_name_lst"),
              "Select Data Category",
              choices = c("variants", "annotations", "genotypes"),
              selected = "genotypes",
              width = "50%"
            ),
            reactable::reactableOutput(ns("results_lst"))
          )
        )
      )
    )
  }

 # ------------------------------------------------------------------
# Cmbined Coordinate Entry Card
# ------------------------------------------------------------------
coordinate_entry_card <- function(ns) {
  bslib::navset_card_underline(
    id = ns("coord_entry_tabs"),
    title = tags$span(
      class = "fw-bold",
     # icon("location-crosshairs", class = "me-2 text-primary"),
      "Set Genomic Coordinates: "
    ),

    # ── GFF3 IMPORT ──
    bslib::nav_panel(
      "Import GFF3",
      icon = icon("file-code"),
      div(
        class = "p-2",
        textInput(
          ns("gene_name"),
          tags$span(
            class = "text-muted small text-uppercase fw-bold",
            "Gene Name (Sobic ID)"
          ),
          value = "",
          placeholder = "e.g. Sobic.005G213600",
          width = "100%"
        ),

        # Inline radio buttons save vertical space
        radioButtons(
          ns("input_method"),
          tags$span(
            class = "text-muted small text-uppercase fw-bold",
            "GFF File Source"
          ),
          choices = c("Remote URL" = "url", "Local Upload" = "file"),
          selected = "url",
          inline = TRUE
        ),

        conditionalPanel(
          condition = paste0("input['", ns("input_method"), "'] === 'url'"),
          textInput(
            ns("gff_url"),
            NULL, # Label omitted because the radio button provides context
            width = "100%",
            value = "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
          )
        ),

        conditionalPanel(
          condition = paste0("input['", ns("input_method"), "'] === 'file'"),
          fileInput(
            ns("gff_file"),
            NULL,
            accept = c(
              ".gff3",
              ".gff",
              ".gff3.gz",
              ".gff.gz",
              "application/gzip",
              "application/x-gzip"
            ),
            width = "100%"
          )
        ),

        div(
          class = "d-grid mt-3",
          actionButton(
            ns("submit"),
            "Get Coordinates",
            class = "btn-success fw-bold py-2",
            icon = icon("search")
          )
        )
      )
    ),

    # ── TAB 2: MANUAL ENTRY ──
    bslib::nav_panel(
      "Manual Entry",
      icon = icon("keyboard"),
      div(
        class = "p-2",
        numericInput(
          ns("chrom"),
          tags$span(
            class = "text-muted small text-uppercase fw-bold",
            "Chromosome Number"
          ),
          value = NULL,
          min = 1,
          step = 1,
          width = "100%"
        ),

        bslib::layout_columns(
          col_widths = c(6, 6),
          numericInput(
            ns("start"),
            tags$span(
              class = "text-muted small text-uppercase fw-bold",
              "Start Position"
            ),
            value = NULL,
            min = 1,
            step = 1,
            width = "100%"
          ),
          numericInput(
            ns("end"),
            tags$span(
              class = "text-muted small text-uppercase fw-bold",
              "End Position"
            ),
            value = NULL,
            min = 1,
            step = 1,
            width = "100%"
          )
        ),

        div(
          class = "d-grid mt-3",
          actionButton(
            ns("set_genocod_btn"),
            "Submit Coordinates",
            class = "btn-success fw-bold py-2",
            icon = icon("check")
          )
        )
      )
    )
  )
}

# ------------------------------------------------------------------
# Main Gene Coordinates Tab
# ------------------------------------------------------------------
gene_cord_tab <- function(ns) {
  tagList(
    uiOutput(ns("genomic_range_vboxes")),

    tags$hr(class = "my-4"),

    # Use a navset to hold both the hotspot plot and the new LD analysis
    bslib::navset_card_tab(
      id = ns("analysis_tabs"),
      title = tags$span(
        class = "fw-bold",
        "Decision Support: "
      ),
      full_screen = TRUE,
      bslib::nav_panel(
        "Variant Hotspots",
        icon = icon("fire"),
        uiOutput(outputId = ns("variant_hotspot_plot_ui"))
      )
    ),

    # Gff loading card
    bslib::layout_columns(
      col_widths = c(-3, 6, -3),
      coordinate_entry_card(ns)
    )
  )
}

  query_actions_tab <- function(ns) {
    tagList(
      bslib::layout_columns(
        col_widths = c(-3,6,-3),
        query_database_card(ns)
      ),
      bslib::navset_card_tab(
        id = ns("query_db_nav_id"),
        selected = "Main Database Results",
        bslib::nav_panel(
          "Main Database Results",
          uiOutput(ns("query_db_display"))
        )
      )
    )
  }

  get_pcv_card <- function(ns) {
    tagList(
      bslib::layout_columns(
        col_widths = c(-2, 8, -2),
        query_action_card(ns)
      ),
      bslib::navset_card_tab(
        id = ns("pcv_nav_id"),
        selected = "PCVs for KASP Marker Design",
          bslib::nav_panel(
            "PCVs for KASP Marker Design",

            # Accordion for Metadata Filtering
            bslib::accordion(
              id = ns("metadata_filter_accordion"),
              open = FALSE, 
              bslib::accordion_panel(
                "Filter by Sample Metadata",
                icon = icon("filter"),
                # Allow selecting multiple metadata columns
                selectizeInput(
                  ns("meta_filter_col"),
                  "Select Metadata Columns to Filter By:",
                  choices = NULL,
                  multiple = TRUE,
                  width = "100%"
                ),
                # This UI will be dynamically generated
                uiOutput(ns("meta_filter_values_ui")),
                
                # Put start column and buttons on the same row
                fluidRow(
                  column(6, 
                         numericInput(ns("genotype_start_col_input"), "Genotype Start Column:", value = 11, min = 1, step = 1, width = "100%")
                  ),
                  column(3, 
                         actionButton(ns("apply_meta_filter"), "Apply Filters", icon = icon("check"), class = "btn-primary w-100", style = "margin-top: 25px;")
                  ),
                  column(3, 
                         actionButton(ns("clear_meta_filter"), "Clear Filters", icon = icon("xmark"), class = "btn-outline-secondary w-100", style = "margin-top: 25px;")
                  )
                ),
                tags$hr(),
                tags$small(textOutput(ns("filtered_sample_count")), class = "text-muted")
              ),
              class = "mb-3" # Add some margin below the accordion
            ),
            uiOutput(ns("pcvs_kasp_marker_design_result"))
          ),

          bslib::nav_panel(
            "Linkage Disequilibrium",
            icon = icon("link"),
            uiOutput(ns("ld_analysis_ui"))
          )
      )
    )
  }

  query_database_card <- function(ns) {
    bslib::card(
      class = "shadow p",
      bslib::card_header(
        tags$strong("Query Genomic Data for Target Region"),
        class = 'text-center',
        style = "font-size:18px;"
      ),
      selectInput(
        width = "100%",
        inputId = ns("query_database"),
        label = "Select genomic view to explore:",
        choices = c(
          "Genotypes" = "genotypes",
          "Variants" = "variants",
          "Annotations" = "annotations",
          "Annotation Summary" = "annotation_summary"
        ),
        selected = "annotations"
      ),
      # Conditionally show gene name input only when 'Annotations' is selected
      conditionalPanel(
        condition = paste0("input['", ns("query_database"), "'] === 'annotations'"),
        textInput(
          ns("query_gene_name"),
          "Gene Name (Optional)",
          value = NULL,
          placeholder = "e.g. Sobic.005G213600",
          width = "100%"
        )
      ),
      bslib::card_footer(
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            ns("query_dbase_btn"),
            tags$b("Fetch Data"),
            width = "70%",
            icon = icon("play"),
            style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
            `onmouseover` = "this.style.backgroundColor='#145214'",
            `onmouseout` = "this.style.backgroundColor='forestgreen'"
          )
        )
      )
    )
  }

  query_action_card <- function(ns) {
    bslib::card(
      class = "shadow p",
      bslib::card_header(
        tags$strong("Extract Putative Causal Variants"),
        class = 'text-center',
        style = "font-size:18px;"
      ),
      bslib::navset_card_underline(
        id = ns('impact_card'),
        bslib::nav_panel(
          "By Impact & AF",
          icon = icon("filter"),
          bslib::card_body(
            uiOutput(ns("impact_level_ui")),
            sliderInput(
              ns("af_range"),
              "Minimum Alternate Allele Frequency",
              min = 0,
              max = 1,
              value = 0.05,
              step = 0.01,
              width = "100%"
            )
           # verbatimTextOutput(outputId = ns("alt_freq_range"))
          )
        ),
        bslib::nav_panel(
          title = "By Variant ID",
          icon = icon("list-check"),
          bslib::card_body(
            selectizeInput(
              ns("manual_variant_ids"),
              "Select Variant IDs",
              choices = NULL,
              multiple = TRUE,
              width = "100%",
              options = list(placeholder = "Select variants from region...")
            )
          )
        )
      ),
      bslib::card_footer(
        div(
          class = "d-grid",
          actionButton(
            ns("get_pcv_btn"),
            tags$b("Extract PCVs"),
            class = "btn-warning",
            width = "100%",
            icon = icon("filter")
          )
        )
      )
    )
  }

  # ------------------------------------------------------------------
  # THE TWO MAIN VIEWS THAT'S PRE-CONNECTION AND ACTIVE-CONNECTION UI
  # ------------------------------------------------------------------

  # PRE-CONNECTION  UI
  pre_connection_view <- div(
    class = "container py-5",
    tags$h3(
      "Choose database connection mode",
      class = "text-center fw-bold mb-4"
    ),
    bslib::layout_columns(
      col_widths = c(2, 4, 4, 2),
      div(),
      bslib::card(
        class = "text-center shadow-sm p-4",
        style = "border-top: 4px solid #3498DB;",
        tags$div(class = "h1 text-info mb-3", icon("database")),
        tags$h5("Local Database Folder (Parquet)", class = "fw-bold"),
        p(
          "Perform offline variant discovery using a local Parquet directory.",
          class = "text-muted mb-4"
        ),
        shinyFiles::shinyDirButton(
          id = ns("connect_btn"),
          label = "Browse for Parquet Folder",
          icon = icon("folder-open"),
          title = "Select Database Folder",
          class = "btn-info btn-lg rounded-pill px-4",
          multiple = FALSE
        )
      ),
      bslib::card(
        class = "text-center shadow-sm p-4",
        style = "border-top: 4px solid #27AE60;",
        tags$div(class = "h1 text-success mb-3", icon("server")),
        tags$h5("Cloud / Remote Database", class = "fw-bold"),
        p(
          "Perform online variant discovery using a remote curated sorghum variant resource.",
          class = "text-muted mb-4"
        ),
        actionButton(
          ns("btn_show_postgres_modal"),
          "Connect to Server",
          class = "btn-success btn-lg rounded-pill px-4",
          icon = icon("plug")
        )
      ),
      div()
    )
  )
  

  active_dashboard_view <- bslib::layout_sidebar(
    fillable = FALSE,
    sidebar = bslib::sidebar(
      id = ns('db_sidebar'),
      width = 320,
      bg = "#f8f9fa",
      class = "p-3",
      title = tags$h3(
        "MAIN MENU",
        class = "text-muted fw-bold mb-3",
        style = "letter-spacing: 1.5px; font-size: 0.9rem;"
      ),

      div(
        class = "d-flex align-items-center justify-content-center px-3 py-2 mb-4 bg-success-subtle text-success rounded-pill border border-success-subtle shadow-sm",
        icon("wifi", class = "me-2"),
        tags$span("Database Connected", class = "fw-bold small mb-0")
      ),

      # Navigation Buttons
      # Get database Info button
      actionButton(
        ns("get_db_info"),
        "Database Overview",
        icon = icon("info-circle"),
        width = "100%",
        class = "btn-primary btn-lg mb-2 fw-bold text-start"
      ),

      # Set Gene Cordinates Button
      actionButton(
        ns("show_gene_cord_btn"),
        "Target Region",
        icon = icon("crosshairs"),
        width = "100%",
        class = "btn-outline-primary btn-lg mb-4 fw-bold text-start"
      ),

      # Query Database button
      actionButton(
        ns("show_query_actions_btn"),
        "Browse Variants",
        icon = icon("search"),
        width = "100%",
        class = "btn-outline-primary btn-lg mb-4 fw-bold text-start"
      ),

      # Get Putative Causal Variants button
      actionButton(
        ns("get_pcv_sidebar_btn"),
        "Causal Variants",
        icon = icon("bolt"),
        width = "100%",
        class = "btn-outline-primary btn-lg mb-4 fw-bold text-start"
      ),

      # Design KASP Markers button
      actionButton(
        ns("design_kasp_sidebar_btn"),
        "Marker Design",
        icon = icon("dna"),
        width = "100%",
        class = "btn-outline-primary btn-lg mb-4 fw-bold text-start"
      ),

      actionButton(
        ns("disconnect_btn"),
        "Disconnect",
        icon = icon("power-off"),
        class = "btn-danger w-100 mb-4 py-2 fw-bold"
      )
    ),

    # MAIN CONTENT AREA
    bslib::navset_hidden(
      id = ns('param_header'),

      # Info tab
      bslib::nav_panel_hidden(
        value = 'info_tab',
        info_tab(ns)
      ),

      # Gene coordinates tab
      bslib::nav_panel_hidden(
        value = 'gene_cord',
        gene_cord_tab(ns)
      ),

      # Query actions tab
      bslib::nav_panel_hidden(
        value = 'query_tab',
        query_actions_tab(ns)
      ),

      # Filter out Pcvs
      bslib::nav_panel_hidden(
        value = 'pcv_tab',
        get_pcv_card(ns)
      ),

      # HIDDEN KASP MARKER DESIGN TAB
      bslib::nav_panel_hidden(
        value = "mark_design",

        # Page header with a modern back button
        div(
          class = "d-flex align-items-center mb-4",
          actionButton(
            ns("go_back"),
            "Back to Causal Variants",
            icon = icon("arrow-left"),
            class = "btn btn-light btn-lg me-3" # A more modern, subtle button
          ),
          tags$h4("KASP Marker Design", class = "m-0 fw-bold text-primary")
        ),

        fluidRow(
          column(
            width = 4,
            bslib::card(
              class = "shadow-sm border-0 mb-3",
              bslib::card_body(
                fileInput(
                  ns("modal_genome_file"),
                  "Genome Reference File",
                  accept = c(".fa", ".fasta", ".gz"),
                  width = "100%"
                ),
                selectizeInput(
                  ns("modal_marker_ID"),
                  "Marker ID",
                  choices = NULL,
                  options = list(placeholder = "Select variants..."),
                  multiple = TRUE,
                  width = "100%"
                ),
                textInput(
                  ns("modal_reg_name"),
                  "Region Name",
                  width = "100%",
                  placeholder = "lgs1"
                ),
                numericInput(
                  ns("modal_maf"),
                  "Minor Allele Frequency (MAF)",
                  value = 0.05,
                  min = 0,
                  max = 1,
                  step = 0.01,
                  width = "100%"
                ),
                numericInput(
                  ns("modal_genotype_start_col"),
                  "Index of First Genotype Column",
                  value = 11,
                  min = 1,
                  step = 1,
                  width = "100%"
                ) 
              ),
              bslib::card_footer(
                class = "bg-transparent border-0",
                div(
                  class = "mt-2 d-grid gap-2",
                  actionButton(
                    ns("modal_run_but"),
                    "Design KASP Marker",
                    icon = icon("play", class = "me-2"),
                    class = "btn-success btn-lg py-3 fw-bold shadow-sm"
                  )
                )
              )
            )
          ),

          # Right Column: Results Tabs
          column(
            width = 8,
            bslib::navset_card_pill(
              id = ns("results_tabs"),
              full_screen = TRUE,

              # Marker Data Sub-tab
              bslib::nav_panel(
                title = "Marker Data",
                icon = icon("table"),
                value = "table_tab",
                bslib::card(
                  class = "shadow-sm border-0",
                  bslib::card_header(
                    class = "bg-light d-flex align-items-center",
                    icon("table", class = "me-2 text-primary"),
                    strong("KASP Marker Design Results & Sequence Alignment")
                  ),
                  bslib::card_body(
                    class = "p-0",
                    DT::DTOutput(ns("kasp_table"), height = "500px")
                  ),
                  bslib::card_footer(
                    class = "bg-light",
                    bslib::layout_columns(
                      col_widths = c(3, 6, 3),
                      selectInput(
                        ns("exten"),
                        "Download file as?",
                        choices = c(".csv", ".xlsx"),
                        selected = ".xlsx"
                      ),
                      textInput(
                        ns("file_name"),
                        "Enter File Prefix",
                        value = "Kasp M_D for Intertek"
                      ),
                      div(
                        class = "d-flex align-items-end h-100 pb-3",
                        downloadButton(
                          ns("download_table"),
                          "Export",
                          class = "btn-success w-100",
                          icon = icon("download")
                        )
                      )
                    )
                  )
                )
              ),

              # Alignment Plot Sub-tab
              bslib::nav_panel(
                title = "Alignment Plot",
                icon = icon("chart-bar"),
                value = "plot_tab",
                bslib::card(
                  class = "shadow-sm border-0",
                  bslib::card_header(
                    class = "bg-light d-flex align-items-center",
                    icon("chart-line", class = "me-2 text-primary"),
                    strong("Sequence Alignment Visualization")
                  ),
                  bslib::card_body(
                    tagList(
                      selectizeInput(
                        ns("plot_choice"),
                        "Select Marker ID",
                        width = "30%",
                        choices = NULL
                      ),
                      uiOutput(ns("plot_container")),
                      plotOutput(ns("plot"), height = "400px"),
                      downloadButton(
                        ns("download_plot"),
                        "Export All Plots (PDF)",
                        class = "btn-success mt-4 px-4",
                        icon = icon("download")
                      )
                    )
                  ),
                  bslib::card_footer(
                    class = "bg-light text-muted small",
                    icon("info-circle", class = "me-1"),
                    "Interactive alignment plot showing variant positions relative to reference genome"
                  )
                )
              )
            )
          )
        )
      )
    )
  )



  # ------------------------------------------------------------------
  #  MAIN UI TAGLIST 
  # ------------------------------------------------------------------
  tagList(
    shinyjs::useShinyjs(),
    div(id = ns("pre_connection_panel"), pre_connection_view),
    shinyjs::hidden(div(
      id = ns("active_dashboard_panel"),
      active_dashboard_view
    ))
  )
}




#' Variant Discovery Server Function
#'
#' @description Server logic for variant discovery module.
#'
#' @param id Internal parameter for {shiny}.
#' @importFrom shiny showModal modalDialog
#'
#' @noRd
mod_variant_discovery_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ------------------------------------------------------------------
    # DASHBOARD  SIDEBAR
    # ------------------------------------------------------------------

    # database info
    observeEvent(input$get_db_info, {
      updateTabsetPanel(session, "param_header", selected = "info_tab")
      update_sidebar_buttons("get_db_info")
    })

    # Gene cordinates
    observeEvent(input$show_gene_cord_btn, {
      updateTabsetPanel(session, "param_header", selected = "gene_cord")
      update_sidebar_buttons("show_gene_cord_btn")
    })

    # Query database
    observeEvent(input$show_query_actions_btn, {
      updateTabsetPanel(session, "param_header", selected = "query_tab")
      update_sidebar_buttons("show_query_actions_btn")
    })

    # Design KASP Markers (from sidebar)
    observeEvent(input$design_kasp_sidebar_btn, {
      # Evaluate target genotype matrix to isolate metadata-filtered accessions prior to UI rendering
      active_data <- if (!is.null(values$filtered_pcvs_by_meta)) {
        values$filtered_pcvs_by_meta
      } else {
        values$query_geno_react
      }
      req(active_data)

      bslib::toggle_sidebar(id = 'db_sidebar', open = 'closed')
      updateTabsetPanel(inputId = 'param_header', selected = "mark_design")
      updateSelectizeInput(
        session,
        inputId = "modal_marker_ID",
        choices = active_data$variant_id,
        server = TRUE
      )
      update_sidebar_buttons("design_kasp_sidebar_btn")
    })

    # ------------------------------------------------------------------
    # REACTIVE VALUES & STATE MANAGEMENT
    # ------------------------------------------------------------------
    rv <- reactiveValues(
      conn = NULL,
      connected = FALSE,
      conn_type = NULL,
      tables = NULL,
      status_message = "Not connected",
      db_path = NULL,
      sample_metadata = NULL,
      # Data stores
      variant_impact = NULL,
      sqlite_summary = NULL,
      variant_count = NULL,
      variant_stats = NULL,
      lst_tbl_column = NULL
    )

    output$is_connected <- reactive(rv$connected)
    outputOptions(output, "is_connected", suspendWhenHidden = FALSE)

    # ------------------------------------------------------------------
    # DATABASE CONNECTION MODE LOGIC
    # ------------------------------------------------------------------

    #  Local DuckDB/Parquet Connection
    volumes <- c(Home = path.expand("~"), shinyFiles::getVolumes()())
    shinyFiles::shinyDirChoose(
      input,
      "connect_btn",
      roots = volumes,
      session = session
    )

    observeEvent(input$connect_btn, {
      shinyjs::hide(id = 'impact_card')

      # Check if a valid directory was selected
      if (!is.integer(input$connect_btn)) {
        db_path <- shinyFiles::parseDirPath(volumes, input$connect_btn)

        shinybusy::show_modal_spinner(
          spin = "fading-circle",
          color = "#27AE60",
          text = "Connecting to Database... Please wait."
        )
        tryCatch(
          {
            if (!is.null(rv$conn) && rv$conn_type == "sqlite") {
              disconnect_local_db(rv$conn)
              rv$conn <- NULL
            }

            rv$conn <- connect_local_db(folder_path = db_path)
            rv$connected <- TRUE
            rv$conn_type <- "sqlite"
            rv$tables <- list_sqlite_tables(con = rv$conn)
            rv$sample_metadata <- get_sample_metadata(con = rv$conn)
            rv$db_path <- db_path

            shinyjs::hide("pre_connection_panel")
            shinyjs::show("active_dashboard_panel")
            shinyWidgets::show_alert(
              title = "Success!",
              text = "Local Database connected successfully",
              type = "success",
              timer = 3000
            )
          },
          error = function(e) {
            shinyWidgets::show_alert(
              title = "Connection Error!",
              text = paste(
                "Failed to connect or read from database:",
                e$message
              ),
              type = "danger",
              timer = 8000
            )
          },
          finally = {
            shinybusy::remove_modal_spinner()
          }
        )
      }
    })

    # Show Modal for PostgreSQL Connection
    observeEvent(input$btn_show_postgres_modal, {
      showModal(
        modalDialog(
          title = div(
            class = "d-flex align-items-center",
            icon(
              "cloud",
              class = "text-success me-2",
              style = "font-size: 1.5rem;"
            ),
            tags$b("Connect to Cloud API")
          ),
          size = "m",
          easyClose = TRUE,
          fade = TRUE,
          div(
            class = "mb-4 text-muted",
            "By default, panGB connects to the public ",
            tags$b("Curated Sorghum Variant Resource"),
            " database. ",
            "If your institution hosts a custom database for another crop, you can specify the endpoint below."
          ),
          radioButtons(
            ns("api_choice"),
            label = tags$b("Select Database Source:"),
            choices = c(
              "Curated Sorghum Variant Resource" = "default",
              "Custom API Endpoint" = "custom"
            ),
            selected = "default"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'custom'", ns("api_choice")),
            textInput(
              ns("api_url"),
              label = "Server Address (URL)",
              placeholder = "e.g., http://rice-server:8000",
              width = "100%"
            )
          ),
          footer = tagList(
            modalButton("Cancel"),
            actionButton(
              ns("btn_connect_postgres"),
              "Connect",
              class = "btn-success fw-bold px-4",
              icon = icon("plug")
            )
          )
        )
      )
    })

    observeEvent(input$btn_connect_postgres, {
      removeModal()
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Connecting to database..."
      )
      tryCatch(
        {
          target_url <- if (input$api_choice == "custom") {
            trimws(input$api_url)
          } else {
            get_api_url()
          }

          set_api_url(target_url)
          rv$tables <- pg_list_tables()
          rv$connected <- TRUE
          rv$sample_metadata <- pg_get_sample_metadata()
          rv$conn_type <- "postgres"

          shinyjs::hide("pre_connection_panel")
          shinyjs::show("active_dashboard_panel")
          shinyWidgets::show_alert(
            title = "Success!",
            text = "Database connected successfully",
            type = "success",
            timer = 3000
          )
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Failed!",
            text = "Unable to reach database. Please check your internet connection",
            type = "danger",
            timer = 5000
          )
        },
        finally = {
          shinybusy::remove_modal_spinner()
        }
      )
    })

    # Disconnect
    observeEvent(input$disconnect_btn, {
      if (!is.null(rv$conn) && rv$conn_type == "sqlite") {
        disconnect_local_db(rv$conn)
      }
      rv$conn <- NULL
      rv$connected <- FALSE
      rv$conn_type <- NULL
      rv$tables <- NULL
      rv$db_path <- NULL
      rv$sample_metadata <- NULL

      rv$variant_impact <- NULL
      rv$sqlite_summary <- NULL
      rv$variant_count <- NULL
      rv$variant_stats <- NULL
      rv$lst_tbl_column <- NULL

      values$result <- NULL
      values$query_db_val <- NULL
      values$query_ann_react <- NULL
      values$query_geno_react <- NULL

      shinyjs::show("pre_connection_panel")
      shinyjs::hide("active_dashboard_panel")
      shinyWidgets::show_toast(
        title = "Disconnected from Database",
        type = "warning",
        timer = 3000
      )
    })

    session$onSessionEnded(function() {
      isolate({
        if (!is.null(rv$conn) && rv$conn_type == "sqlite") {
          disconnect_local_db(rv$conn)
          rv$conn <- NULL
        }
      })
    })

    # ------------------------------------------------------------------
    #   METADATA MAP LOGIC
    # ------------------------------------------------------------------
    observeEvent(rv$sample_metadata, {
      req(rv$sample_metadata)

      meta <- rv$sample_metadata
      all_cols <- colnames(meta)

      preferred_cols <- c(
        "countryorigin",
        "population",
        "region",
        "germplasm",
        "botanicaltype",
        "kmeans_cluster"
      )

      cols_to_exclude <- c(
        "lat",
        "lon",
        "latitude",
        "lib",
        "sample",
        "plantname",
        "pinumber",
        "array_index",
        "accession_id",
        "sample_id",
        "isnum",
        "wasap",
        "duplicatelibdecision",
        "filtermissing"
      )

      actual_preferred <- intersect(preferred_cols, all_cols)
      candidate_cols <- setdiff(all_cols, c(actual_preferred, cols_to_exclude))

      heuristic_passed_cols <- Filter(
        function(col_name) {
          col_data <- meta[[col_name]]

          if (all(is.na(col_data))) {
            return(FALSE)
          }

          valid_data <- na.omit(col_data)
          if (length(valid_data) == 0) {
            return(FALSE)
          }

          num_unique <- length(unique(valid_data))

          if (num_unique <= 1) {
            return(FALSE)
          }

          if (num_unique > 0.9 * length(valid_data) && num_unique > 100) {
            return(FALSE)
          }

          if (is.numeric(col_data)) {
            return(num_unique < 25)
          }

          if (is.character(col_data) || is.factor(col_data)) {
            return(num_unique < 100)
          }

          return(FALSE)
        },
        candidate_cols
      )

      suitable_cols <- c(actual_preferred, heuristic_passed_cols)

      updateSelectInput(
        session,
        "color_by_col",
        choices = suitable_cols,
        selected = if ("countryorigin" %in% suitable_cols) {
          "countryorigin"
        } else {
          suitable_cols[1]
        }
      )
    })

    output$accession_map <- leaflet::renderLeaflet({
      req(rv$sample_metadata, input$color_by_col, rv$conn_type)

      if (rv$conn_type == "sqlite") {
        query_map_accessions(
          metadata = rv$sample_metadata,
          color_by = input$color_by_col
        )
      } else {
        pg_map_accessions(
          metadata = rv$sample_metadata,
          color_by = input$color_by_col
        )
      }
    })

    # ------------------------------------------------------------------
    # DATABASE INFO TAB
    # ------------------------------------------------------------------

    fetch_data <- function(
      button_id,
      task_name,
      data_slot,
      sqlite_func_name,
      pg_func_name
    ) {
      observeEvent(input[[button_id]], {
        req(rv$connected)

        shinyWidgets::show_toast(
          title = "Fetching Data",
          text = paste("Retrieving", task_name, "in the background..."),
          type = "info",
          timer = 3000,
          position = "bottom-end"
        )

        c_type <- rv$conn_type
        d_path <- rv$db_path

        p <- future::future(
          {
            if (c_type == "sqlite") {
              temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
              res <- do.call(sqlite_func_name, list(con = temp_con))
              disconnect_local_db(temp_con, quiet = TRUE)
              res
            } else {
              do.call(pg_func_name, list())
            }
          },
          seed = TRUE,
          packages = c("panGenomeBreedr")
        )

        promises::then(
          p,
          onFulfilled = function(res) {
            rv[[data_slot]] <- res
            shinyWidgets::show_toast(
              title = "Complete",
              text = paste(task_name, "updated."),
              type = "success",
              timer = 3000,
              position = "bottom-end"
            )
          },
          onRejected = function(err) {
            shinyWidgets::show_alert(
              title = "Failed!",
              text = err$message,
              type = "danger"
            )
          }
        )
      })
    }

    fetch_data(
      "get_impact",
      "Variant Impact Summary",
      "variant_impact",
      "variant_impact_summary",
      "pg_variant_impact_summary"
    )

    fetch_data(
      "get_stats",
      "Variant Statistics",
      "variant_stats",
      "variant_stats",
      "pg_variant_stats"
    )

    fetch_data(
      "get_summary",
      "Table Summaries",
      "sqlite_summary",
      "summarize_sqlite_tables",
      "pg_summarize_tables"
    )

    fetch_data(
      "get_types",
      "Variant Type Counts",
      "variant_count",
      "count_variant_types",
      "pg_count_variant_types"
    )

    observeEvent(input$get_schema, {
      req(rv$connected, input$table_name_lst)

      shinyWidgets::show_toast(
        title = "Fetching Data",
        text = paste(
          "Retrieving schema for",
          input$table_name_lst,
          "in the background..."
        ),
        type = "info",
        timer = 3000,
        position = "bottom-end"
      )

      c_type <- rv$conn_type
      d_path <- rv$db_path
      t_name <- input$table_name_lst

      p <- future::future(
        {
          if (c_type == "sqlite") {
            temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
            res <- list_table_columns(con = temp_con, table_name = t_name)
            disconnect_local_db(temp_con, quiet = TRUE)
            res
          } else {
            pg_list_table_columns(table_name = t_name)
          }
        },
        seed = TRUE,
        packages = c("panGenomeBreedr")
      )

      promises::then(
        p,
        onFulfilled = function(res) {
          rv$lst_tbl_column <- res
          shinyWidgets::show_toast(
            title = "Complete",
            text = "Table schema updated.",
            type = "success",
            timer = 3000,
            position = "bottom-end"
          )
        },
        onRejected = function(err) {
          shinyWidgets::show_alert(
            title = "Failed!",
            text = err$message,
            type = "danger"
          )
        }
      )
    })

    output$table_impact_id <- reactable::renderReactable({
      req(rv$variant_impact)
      render_reactable(rv$variant_impact)
    })
    output$table_var_stats_id <- reactable::renderReactable({
      req(rv$variant_stats)
      render_reactable(rv$variant_stats)
    })
    output$sum_sqlite_id <- reactable::renderReactable({
      req(rv$sqlite_summary)
      render_reactable(rv$sqlite_summary)
    })
    output$count_variant_typ_id <- reactable::renderReactable({
      req(rv$variant_count)
      render_reactable(rv$variant_count)
    })
    output$results_lst <- reactable::renderReactable({
      req(rv$lst_tbl_column)
      render_reactable(rv$lst_tbl_column)
    })

    # ------------------------------------------------------------------
    # QUERY LOGIC
    # ------------------------------------------------------------------
    values <- reactiveValues(
      result = NULL,
      query_db_val = NULL,
      query_ann_react = NULL,
      query_geno_react = NULL,
      last_action = NULL,
      regional_variants = NULL,
      regional_annotations = NULL,
      regional_genotypes = NULL,
      ld_results = NULL,
      filtered_pcvs_by_meta = NULL,
      geodesic_plot_obj = NULL
    )

    show_toast_success <- function(text, type = "success", timer = 3000) {
      shinyWidgets::show_toast(
        title = "Success",
        text = text,
        timer = timer,
        position = "bottom-end",
        type = type
      )
    }

    annotation_summary_results_ui <- function(ns) {
      bslib::card(
        bslib::card_header("Annotation Statistics"),
        bslib::card(
          bslib::card_header(
            "Annotation Summary",
            class = "bg-primary text-light"
          ),
          reactable::reactableOutput(ns("ann_summary_tbl"))
        ),
        bslib::card(
          bslib::card_header("Impact Summary", class = "bg-primary text-light"),
          reactable::reactableOutput(ns("impact_summary_tbl"))
        ),
        bslib::card(
          bslib::card_header(
            "Variant Type Totals",
            class = "bg-primary text-light"
          ),
          reactable::reactableOutput(ns("variant_totals_tbl"))
        ),
        class = "mt-3"
      )
    }

    observeEvent(values$result, {
      req(rv$connected, values$result)

      values$regional_variants <- NULL
      values$regional_annotations <- NULL
      values$regional_genotypes <- NULL
      values$ld_results <- NULL
      values$geodesic_plot_obj <- NULL

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Plotting variants hotspot..."
      )

      c_type <- rv$conn_type
      d_path <- rv$db_path
      chr <- values$result$chrom
      st <- values$result$start
      en <- values$result$end

      p <- future::future(
        {
          fetch_func <- if (c_type == "sqlite") query_db else pg_query_db
          con_arg <- if (c_type == "sqlite") {
            list(con = connect_local_db(d_path, quiet = TRUE))
          } else {
            list()
          }
          on.exit(
            if (c_type == "sqlite") {
              disconnect_local_db(con_arg$con, quiet = TRUE)
            }
          )

          list(
            variants = do.call(
              fetch_func,
              c(
                con_arg,
                list(table_name = "variants", chrom = chr, start = st, end = en)
              )
            ),
            annotations = do.call(
              fetch_func,
              c(
                con_arg,
                list(
                  table_name = "annotations",
                  chrom = chr,
                  start = st,
                  end = en
                )
              )
            ),
            genotypes = do.call(
              fetch_func,
              c(
                con_arg,
                list(
                  table_name = "genotypes",
                  chrom = chr,
                  start = st,
                  end = en
                )
              )
            )
          )
        },
        seed = TRUE,
        packages = c("panGenomeBreedr")
      )

      promises::then(
        p,
        onFulfilled = function(data) {
          values$regional_variants <- data$variants
          values$regional_annotations <- data$annotations
          values$regional_genotypes <- data$genotypes
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_toast("Variant hotspots plotted.", type = "success")
        },
        onRejected = function(err) {
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_alert(
            title = "Data Fetch Error",
            text = err$message,
            type = "danger"
          )
        }
      )
    })

    # ------------------------------------------------------------------
    # METADATA FILTERING LOGIC FOR PCVs
    # ------------------------------------------------------------------

    observeEvent(values$query_geno_react, {
      req(rv$sample_metadata)
      meta_cols <- colnames(rv$sample_metadata)
      filterable_cols <- setdiff(meta_cols, c("lib", "array_index"))
      updateSelectizeInput(
        session,
        "meta_filter_col",
        choices = filterable_cols,
        selected = ""
      )
    })

    output$meta_filter_values_ui <- renderUI({
      req(input$meta_filter_col)

      lapply(input$meta_filter_col, function(col) {
        unique_values <- sort(unique(na.omit(rv$sample_metadata[[col]])))
        selectizeInput(
          ns(paste0("meta_vals_", col)),
          label = paste("Select values for:", col),
          choices = unique_values,
          multiple = TRUE,
          width = "100%"
        )
      })
    })

    observeEvent(input$apply_meta_filter, {
      req(values$query_geno_react, rv$sample_metadata, input$meta_filter_col)

      filters_list <- lapply(input$meta_filter_col, function(col) {
        input[[paste0("meta_vals_", col)]]
      })
      names(filters_list) <- input$meta_filter_col

      filters_list <- filters_list[sapply(filters_list, function(x) {
        !is.null(x) && length(x) > 0
      })]

      if (length(filters_list) == 0) {
        shinyWidgets::show_alert(
          title = "No Filter Values",
          text = "Please select at least one value to filter by.",
          type = "warning"
        )
        return()
      }

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Applying metadata filters..."
      )

      tryCatch(
        {
          genotype_start_col_val <- input$genotype_start_col_input

          result_df <- if (rv$conn_type == "sqlite") {
            query_geno_by_meta(
              con = rv$conn,
              genotype_matrix = values$query_geno_react,
              genotype_start_col = genotype_start_col_val,
              filters = filters_list
            )
          } else {
            pg_query_geno_by_meta(
              genotype_matrix = values$query_geno_react,
              genotype_start_col = genotype_start_col_val,
              filters = filters_list
            )
          }

          values$filtered_pcvs_by_meta <- result_df

          shinybusy::remove_modal_spinner()
          show_toast_success(text = paste("Filtered PCVs by metadata."))
        },
        error = function(e) {
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_alert(
            title = "Filter Error",
            text = e$message,
            type = "danger"
          )
          values$filtered_pcvs_by_meta <- NULL
        }
      )
    })

    observeEvent(input$clear_meta_filter, {
      updateSelectizeInput(session, "meta_filter_col", selected = "")
      values$filtered_pcvs_by_meta <- NULL
      show_toast_success(text = "Metadata filters cleared.", type = "info")
    })

    output$genotype_results_tbl <- reactable::renderReactable({
      if (!is.null(values$filtered_pcvs_by_meta)) {
        req(values$filtered_pcvs_by_meta)
        render_reactable(values$filtered_pcvs_by_meta)
      } else {
        req(values$query_geno_react)
        render_reactable(values$query_geno_react)
      }
    })

    output$pcvs_kasp_marker_design_result <- renderUI({
      req(values$query_geno_react)
      genotype_results_ui(ns)
    })

    output$filtered_sample_count <- renderText({
      if (!is.null(values$filtered_pcvs_by_meta)) {
        original_samples <- ncol(values$query_geno_react) -
          length(c(
            "variant_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "variant_type",
            "major_allele",
            "minor_allele",
            "major_allele_freq",
            "minor_allele_freq"
          ))
        filtered_samples <- ncol(values$filtered_pcvs_by_meta) -
          length(c(
            "variant_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "variant_type",
            "major_allele",
            "minor_allele",
            "major_allele_freq",
            "minor_allele_freq"
          ))
        paste0(
          "Showing ",
          filtered_samples,
          " samples (out of ",
          original_samples,
          " original samples) after filtering."
        )
      } else if (!is.null(values$query_geno_react)) {
        original_samples <- ncol(values$query_geno_react) -
          length(c(
            "variant_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "variant_type",
            "major_allele",
            "minor_allele",
            "major_allele_freq",
            "minor_allele_freq"
          ))
        paste0(
          "Showing all ",
          original_samples,
          " samples. No metadata filters applied."
        )
      } else {
        "No samples to display."
      }
    })

    observeEvent(input$pcv_nav_id, {
      if (input$pcv_nav_id == "Linkage Disequilibrium") {
        req(values$regional_genotypes)
        #req(values$query_geno_react)

        if (
          nrow(values$regional_genotypes) > 0 &&
            "variant_id" %in% names(values$regional_genotypes)
        ) {
          variant_choices <- values$regional_genotypes$variant_id

            pre_selected_variants <- if (
              !is.null(values$query_geno_react) &&
                nrow(values$query_geno_react) > 0
            ) {
              intersect(values$query_geno_react$variant_id, variant_choices)
            } else {
              NULL
            }

          updateSelectizeInput(
            session,
            "target_ids",
            selected = pre_selected_variants,
            choices = variant_choices,
            server = TRUE
          )
        } else {
          updateSelectizeInput(
            session,
            "target_ids",
            choices = character(0),
            server = TRUE
          )
        }
      }
    })

    observeEvent(input$submit, {
      removeModal()
      req(input$gene_name)
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Getting gene coordinates... Please wait."
      )
      tryCatch(
        {
          if (input$input_method == "url") {
            req(input$gff_url)
            gff_path <- input$gff_url
          } else if (input$input_method == "file") {
            req(input$gff_file)
            gff_path <- input$gff_file$datapath
          }
          values$result <- gene_coord_gff(trimws(input$gene_name), gff_path)
          shinyWidgets::show_alert(
            title = "Found Gene Co-ordinates",
            text = sprintf(
              "Chromosome: %s | Start: %d | End: %d",
              values$result$chrom,
              values$result$start,
              values$result$end
            ),
            type = "success",
            timer = 5000
          )
          shinyjs::show(id = 'impact_card')
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Failed!",
            text = e$message,
            type = "danger",
            timer = 5000
          )
        },
        finally = {
          shinybusy::remove_modal_spinner()
        }
      )
    })

    observeEvent(input$set_genocod_btn, {
      req(input$chrom, input$start, input$end, rv$connected)

      chr_val <- sprintf("Chr%02d", input$chrom)
      st_val <- input$start
      en_val <- input$end

      if (st_val >= en_val) {
        shinyWidgets::show_alert(
          title = "Invalid Range",
          text = "Start position must be less than end position.",
          type = "warning"
        )
        return()
      }

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Validating coordinates... Please wait."
      )

      c_type <- rv$conn_type
      d_path <- rv$db_path

      if (!is.null(rv$variant_stats)) {
        v_stats <- rv$variant_stats
        chr_stats <- v_stats[v_stats$chrom == chr_val, ]

        shinybusy::remove_modal_spinner()

        if (nrow(chr_stats) == 0) {
          shinyWidgets::show_alert(
            title = "Invalid Chromosome",
            text = paste("Chromosome", chr_val, "not found in the database."),
            type = "warning"
          )
          return()
        }

        if (en_val < chr_stats$min_pos || st_val > chr_stats$max_pos) {
          shinyWidgets::show_alert(
            title = "Coordinates Out of Bounds",
            text = sprintf(
              "The specified range is outside the available data for %s (Min: %s, Max: %s).",
              chr_val,
              format(chr_stats$min_pos, big.mark = ","),
              format(chr_stats$max_pos, big.mark = ",")
            ),
            type = "warning"
          )
          return()
        }

        values$result <- list(
          chrom = chr_val,
          start = st_val,
          end = en_val
        )
        shinyWidgets::show_alert(
          title = "Gene Coordinates Set",
          text = sprintf(
            "Chromosome: %s | Start: %d | End: %d",
            chr_val,
            st_val,
            en_val
          ),
          type = "success",
          timer = 5000
        )
        shinyjs::show(id = 'impact_card')
      } else {
        p <- future::future(
          {
            if (c_type == "sqlite") {
              temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
              res <- variant_stats(con = temp_con, include_annotations = FALSE)
              disconnect_local_db(temp_con, quiet = TRUE)
              res
            } else {
              pg_variant_stats(include_annotations = FALSE)
            }
          },
          seed = TRUE,
          packages = c("panGenomeBreedr")
        )

        promises::then(
          p,
          onFulfilled = function(v_stats) {
            shinybusy::remove_modal_spinner()

            chr_stats <- v_stats[v_stats$chrom == chr_val, ]

            if (nrow(chr_stats) == 0) {
              shinyWidgets::show_alert(
                title = "Invalid Chromosome",
                text = paste(
                  "Chromosome",
                  chr_val,
                  "not found in the database."
                ),
                type = "warning"
              )
              return()
            }

            if (en_val < chr_stats$min_pos || st_val > chr_stats$max_pos) {
              shinyWidgets::show_alert(
                title = "Coordinates Out of Bounds",
                text = sprintf(
                  "The specified range is outside the available data for %s (Min: %s, Max: %s).",
                  chr_val,
                  format(chr_stats$min_pos, big.mark = ","),
                  format(chr_stats$max_pos, big.mark = ",")
                ),
                type = "warning"
              )
              return()
            }

            values$result <- list(
              chrom = chr_val,
              start = st_val,
              end = en_val
            )
            shinyWidgets::show_alert(
              title = "Gene Coordinates Set",
              text = sprintf(
                "Chromosome: %s | Start: %d | End: %d",
                chr_val,
                st_val,
                en_val
              ),
              type = "success",
              timer = 5000
            )
            shinyjs::show(id = 'impact_card')
          },
          onRejected = function(err) {
            shinybusy::remove_modal_spinner()
            shinyWidgets::show_alert(
              title = "Error Validating Coordinates",
              text = err$message,
              type = "error"
            )
          }
        )
      }
    })

    output$genomic_range_vboxes <- renderUI({
      if (is.null(values$result)) {
        chrom <- "Not Set"
        start_pos <- "Not Set"
        end_pos <- "Not Set"
        region_len <- "Not Set"
      } else {
        chrom <- values$result$chrom
        start_pos <- format(values$result$start, big.mark = ",")
        end_pos <- format(values$result$end, big.mark = ",")

        l <- values$result$end - values$result$start + 1
        if (l < 1000) {
          region_len <- paste0(l, " bp")
        } else if (l < 1e6) {
          region_len <- paste0(round(l / 1000, 2), " kb")
        } else {
          region_len <- paste0(round(l / 1e6, 2), " Mb")
        }
      }

      vbox_title_style <- "font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.5px;"
      vbox_value_style <- "font-size: 1.2rem; font-weight: 700;"
      vbox_icon_style <- "font-size: 1.8rem; opacity: 0.5;"

      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        bslib::value_box(
          title = tags$span("Chromosome", style = vbox_title_style),
          value = tags$span(chrom, style = vbox_value_style),
          showcase = icon("dna", style = vbox_icon_style),
          theme = "white",
          class = "border shadow-sm",
          height = "90px"
        ),
        bslib::value_box(
          title = tags$span("Start Position", style = vbox_title_style),
          value = tags$span(start_pos, style = vbox_value_style),
          showcase = icon("map-marker-alt", style = vbox_icon_style),
          theme = "white",
          class = "border shadow-sm",
          height = "90px"
        ),
        bslib::value_box(
          title = tags$span("End Position", style = vbox_title_style),
          value = tags$span(end_pos, style = vbox_value_style),
          showcase = icon("flag-checkered", style = vbox_icon_style),
          theme = "white",
          class = "border shadow-sm",
          height = "90px"
        ),
        bslib::value_box(
          title = tags$span("Region Length", style = vbox_title_style),
          value = tags$span(region_len, style = vbox_value_style),
          showcase = icon("ruler-horizontal", style = vbox_icon_style),
          theme = "white",
          class = "border shadow-sm",
          height = "90px"
        )
      )
    })

    output$variant_hotspot_plot_ui <- renderUI({
      req(values$result)
      bslib::card(
        class = "shadow-sm border-0",
        bslib::card_header(
          class = "d-flex justify-content-between align-items-center",
          "",
          downloadButton(
            ns("download_hotspot_plot"),
            "Download Plot",
            class = "btn-sm btn-outline-primary"
          )
        ),
        bslib::card_body(
          shinycssloaders::withSpinner(
            plotOutput(ns('hotspot_plot'), height = "600px"),
            type = 4,
            color = "#27AE60"
          )
        )
      )
    })

    hotspot_plot_obj <- reactive({
      req(values$regional_variants, values$result)
      plot_variant_hotspots(
        variant_table = values$regional_variants,
        annotation_table = values$regional_annotations,
        region_start = values$result$start,
        region_end = values$result$end
      )
    })

    output$hotspot_plot <- renderPlot({
      hotspot_plot_obj()
    })

    output$download_hotspot_plot <- downloadHandler(
      filename = function() {
        paste0(
          "hotspot_plot_",
          values$result$chrom,
          "_",
          values$result$start,
          "-",
          values$result$end,
          ".pdf"
        )
      },
      content = function(file) {
        req(hotspot_plot_obj())
        ggplot2::ggsave(
          file,
          plot = hotspot_plot_obj(),
          width = 12,
          height = 8,
          device = "pdf"
        )
      }
    )

    output$ld_analysis_ui <- renderUI({
      req(values$result)
      tagList(
        bslib::card(
          class = "shadow-sm",
          bslib::card_header("LD Analysis Controls"),
          bslib::card_body(
            bslib::layout_columns(
              col_widths = c(3, 3, 3, 3),
              selectInput(
                ns("ld_metric"),
                label = "LD Metric",
                choices = c("R2", "D_prime"),
                selected = "R2"
              ),
              numericInput(
                ns("ld_r2_threshold"),
                label = "Minimum LD Threshold",
                value = 0.2,
                min = 0,
                max = 1,
                step = 0.1
              ),
              numericInput(
                ns("ld_block_threshold"),
                label = "Haploblock Threshold",
                value = 0.8,
                min = 0,
                max = 1,
                step = 0.1
              ),
              selectInput(
                ns("ld_plot_type"),
                label = "Plot Style",
                choices = c("Lines" = FALSE, "Filled Area" = TRUE),
                selected = FALSE
              )
            ),
            bslib::layout_columns(
              col_widths = c(9, 3),
              selectizeInput(
                ns("target_ids"),
                label = "Target Variants (optional)",
                choices = NULL,
                multiple = TRUE,
                width = "100%"
              ),
              selectInput(
                ns("palette_option"),
                label = "Color Palette",
                choices = c(
                  "magma",
                  "plasma",
                  "inferno",
                  "viridis",
                  "mako",
                  "cividis",
                  "rocket",
                  "turbo"
                ),
                selected = "magma"
              )
            )
          ),
          bslib::card_footer(
            actionButton(
              ns("run_ld_analysis"),
              "Run LD Analysis & Plot",
              icon = icon("play"),
              class = "btn-primary w-100"
            )
          )
        ),
        uiOutput(ns("ld_results_output_ui"))
      )
    })

    ld_analysis_results <- eventReactive(input$run_ld_analysis, {
      req(values$regional_genotypes)
      if (nrow(values$regional_genotypes) < 2) {
        shinyWidgets::show_alert(
          title = "Not enough variants",
          text = "LD analysis requires at least two variants in the region.",
          type = "warning"
        )
        return(NULL)
      }
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Computing LD and generating plot..."
      )

      tryCatch(
        {
          target_variants <- if (
            is.null(input$target_ids) || length(input$target_ids) == 0
          ) {
            NULL
          } else {
            input$target_ids
          }

          ld_res <- compute_LD(
            df = values$regional_genotypes,
            target_variant_ids = target_variants
          )

          plot_pkg <- plot_ld_geodesic(
            ld_df = ld_res,
            query_db_geno = values$regional_genotypes,
            query_db_annot = values$regional_annotations,
            metric = input$ld_metric,
            target_variant_ids = target_variants,
            threshold = input$ld_r2_threshold,
            block_threshold = input$ld_block_threshold,
            show_variant_labels = TRUE,
            filled = as.logical(input$ld_plot_type),
            palette_option = input$palette_option
          )

          shinybusy::remove_modal_spinner()
          shinyWidgets::show_toast("LD analysis complete.", type = "success")
          plot_pkg$ld_res <- ld_res
          return(plot_pkg)
        },
        error = function(e) {
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_alert(
            title = "Analysis Error",
            text = e$message,
            type = "error"
          )
          return(NULL)
        }
      )
    })

    output$ld_results_output_ui <- renderUI({
      req(ld_analysis_results())
      bslib::navset_card_pill(
        id = ns("ld_output_tabs"),
        bslib::nav_panel(
          "Geodesic Plot",
          downloadButton(
            ns("download_geodesic_plot"),
            "Download Plot",
            class = "btn-sm btn-outline-primary mb-2"
          ),
          shinycssloaders::withSpinner(plotOutput(
            ns("geodesic_plot_output"),
            height = "600px"
          ))
        ),
        bslib::nav_panel(
          "Haploblock Table",
          downloadButton(
            ns("download_haploblock_table"),
            "Download as Excel",
            class = "btn-sm btn-outline-primary mb-2"
          ),
          shinycssloaders::withSpinner(reactable::reactableOutput(ns(
            "haploblock_table_output"
          )))
        ),
        bslib::nav_panel(
          "LD Results Table",
          downloadButton(
            ns("download_ld_table"),
            "Download as Excel",
            class = "btn-sm btn-outline-primary mb-2"
          ),
          shinycssloaders::withSpinner(reactable::reactableOutput(ns(
            "ld_table_output"
          )))
        )
      )
    })

    output$geodesic_plot_output <- renderPlot({
      results <- ld_analysis_results()
      req(results, results$plot)
      results$plot
    })

    output$haploblock_table_output <- reactable::renderReactable({
      results <- ld_analysis_results()
      req(results)

      tbl <- results$table
      if (!is.null(tbl) && is.data.frame(tbl) && nrow(tbl) > 0) {
        render_reactable(tbl)
      } else {
        reactable::reactable(
          data.frame(Message = character(0)),
          language = reactable::reactableLang(
            noData = "No haploblocks found with the current thresholds."
          )
        )
      }
    })

    output$ld_table_output <- reactable::renderReactable({
      results <- ld_analysis_results()
      req(results, results$ld_res)
      if (
        !is.null(results$ld_res) &&
          is.data.frame(results$ld_res) &&
          nrow(results$ld_res) > 0
      ) {
        render_reactable(results$ld_res)
      } else {
        reactable::reactable(
          data.frame(Message = character(0)),
          language = reactable::reactableLang(
            noData = "No LD results to display."
          )
        )
      }
    })

    output$download_geodesic_plot <- downloadHandler(
      filename = function() {
        paste0(
          "geodesic_plot_",
          values$result$chrom,
          "_",
          values$result$start,
          "-",
          values$result$end,
          ".pdf"
        )
      },
      content = function(file) {
        req(ld_analysis_results()$plot)
        ggplot2::ggsave(
          file,
          plot = ld_analysis_results()$plot,
          width = 14,
          height = 9,
          device = "pdf"
        )
      }
    )

    output$download_haploblock_table <- downloadHandler(
      filename = function() {
        paste0(
          "haploblock_table_",
          values$result$chrom,
          "_",
          values$result$start,
          "-",
          values$result$end,
          ".xlsx"
        )
      },
      content = function(file) {
        if (!requireNamespace("writexl", quietly = TRUE)) {
          shinyWidgets::show_alert(
            title = "Package Missing",
            text = "Please install the 'writexl' package to download Excel files.",
            type = "warning"
          )
          return(NULL)
        }
        req(ld_analysis_results()$table)
        writexl::write_xlsx(ld_analysis_results()$table, file)
      }
    )

    output$download_ld_table <- downloadHandler(
      filename = function() {
        paste0(
          "ld_results_",
          values$result$chrom,
          "_",
          values$result$start,
          "-",
          values$result$end,
          ".xlsx"
        )
      },
      content = function(file) {
        req(ld_analysis_results()$ld_res)
        writexl::write_xlsx(ld_analysis_results()$ld_res, file)
      }
    )

    observe({
      if (
        !is.null(input$query_database) && input$query_database == "q_entire"
      ) {
        output$display_qd_choice <- renderUI({
          tagList(
            selectInput(
              ns("table_name"),
              "Data Type to Extract",
              choices = c("variants", "annotations", "genotypes"),
              selected = "annotations",
              width = "100%"
            ),
            uiOutput(ns("gene_name_id"))
          )
        })
      } else if (
        !is.null(input$query_database) && input$query_database == "q_annt"
      ) {
        output$display_qd_choice <- renderUI({
          tagList(
            selectInput(
              ns("table_name_a"),
              "Annotations Source",
              choices = c("variants", "annotations", "genotypes"),
              selected = "annotations",
              width = "100%"
            ),
            selectInput(
              ns("table_name_v"),
              "Variants Source",
              choices = c("variants", "annotations", "genotypes"),
              selected = "variants",
              width = "100%"
            )
          )
        })
      } else if (is.null(input$query_database)) {
        output$display_qd_choice <- renderUI({
          NULL
        })
      }
    })

    output$gene_name_id <- renderUI({
      if (input$table_name == "annotations") {
        textInput(
          ns("query_gene_name"),
          "Gene Name",
          value = NULL,
          width = "100%"
        )
      }
    })

    # ==========================================================================
    #  QUERY EXECUTION & RESULTS (Browse Variants)
    # ==========================================================================

    observeEvent(input$query_dbase_btn, {
      values$last_action <- NULL
      req(rv$connected, values$result, input$query_database)

      updateTabsetPanel(session, inputId = "query_db_nav_id", selected = "Main Database Results")

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Querying Database... Please wait."
      )

      # Capture reactive inputs for the future
      c_type <- rv$conn_type
      d_path <- rv$db_path
      query_type <- input$query_database
      res_chrom <- values$result$chrom
      res_start <- values$result$start
      res_end <- values$result$end
      q_gene_name <- if (query_type == "annotations" && !is.null(input$query_gene_name) && input$query_gene_name != "") {
        input$query_gene_name
      } else {
        NULL
      }

      p <- future::future({
        # This code runs in a separate process
        if (c_type == "sqlite") {
          temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
          on.exit(disconnect_local_db(temp_con, quiet = TRUE))
        }

        if (query_type %in% c("genotypes", "variants", "annotations")) {
          if (c_type == "sqlite") {
            result <- query_db(
              con = temp_con, table_name = query_type,
              chrom = res_chrom, start = res_start, end = res_end,
              gene_name = q_gene_name
            )
          } else { # postgres
            result <- pg_query_db(
              table_name = query_type,
              chrom = res_chrom, start = res_start, end = res_end,
              gene_name = q_gene_name
            )
          }
          return(list(type = "main_db", data = result))

        } else if (query_type == "annotation_summary") {
          if (c_type == "sqlite") {
            result <- query_ann_summary(
              con = temp_con, variants_table = "variants",
              annotations_table = "annotations", chrom = res_chrom,
              start = res_start, end = res_end
            )
          } else { # postgres
            result <- pg_query_ann_summary(
              annotations_table = "annotations", variants_table = "variants",
              chrom = res_chrom, start = res_start, end = res_end
            )
          }
          return(list(type = "annotation_summ", data = result))
        }
      }, seed = TRUE, packages = c("panGenomeBreedr"))

      promises::then(
        p,
        onFulfilled = function(result) {
          if (result$type == "main_db") {
            values$query_db_val <- result$data
            values$query_ann_react <- NULL
            values$last_action <- "main_db"
            show_toast_success(text = paste("Queried", tools::toTitleCase(query_type)))
          } else if (result$type == "annotation_summ") {
            values$query_ann_react <- result$data
            values$query_db_val <- NULL
            values$last_action <- "annotation_summ"
            show_toast_success("Annotation Summary Retrieved")
          }
          shinybusy::remove_modal_spinner()
        },
        onRejected = function(e) {
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_alert(
            title = "Failed!",
            text = e$message,
            type = "danger",
            timer = 5000
          )
        }
      )
    })

    output$ann_summary_tbl <- reactable::renderReactable({
      req(values$query_ann_react)
      render_reactable(values$query_ann_react$annotation_summary)
    })
    output$impact_summary_tbl <- reactable::renderReactable({
      req(values$query_ann_react)
      render_reactable(values$query_ann_react$impact_summary)
    })
    output$variant_totals_tbl <- reactable::renderReactable({
      req(values$query_ann_react)
      render_reactable(values$query_ann_react$variant_type_totals)
    })

    output$query_db_display <- renderUI({
      req(values$last_action)

      if (values$last_action == "main_db") {
        req(values$query_db_val)
        reactable::renderReactable(render_reactable(values$query_db_val))
      } else if (values$last_action == "annotation_summ") {
        req(values$query_ann_react)
        annotation_summary_results_ui(ns)
      }
    })

    observeEvent(input$get_pcv_sidebar_btn, {
      updateTabsetPanel(session, "param_header", selected = "pcv_tab")
      update_sidebar_buttons("get_pcv_sidebar_btn")
    })

    # ==========================================================================
    # CAUSAL VARIANTS (PCVs) EXTRACTION
    # ==========================================================================
    output$impact_level_ui <- renderUI({
      ns <- session$ns

      impact_order <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
      choices <- character(0)
      selected_choice <- NULL
      placeholder_text <- "First, select a genomic region"

      # Check if regional annotations are available from the hotspot plot data
      if (!is.null(values$regional_annotations) && "impact" %in% names(values$regional_annotations) && nrow(values$regional_annotations) > 0) {
        available_impacts <- unique(na.omit(as.character(values$regional_annotations$impact)))

        if (length(available_impacts) > 0) {
          # Filter and order the available impacts according to our preferred order
          choices <- intersect(impact_order, available_impacts)
          selected_choice <- if (length(choices) > 0) {
            choices[1]
          } else {
            NULL
          } # Default to the highest available impact
          placeholder_text <- "Select impact level(s)..."
        } else {
          placeholder_text <- "No impact data in region"
        }
      }

      # Use selectizeInput to get the placeholder feature
      selectizeInput(
        ns("impact_level"),
        "Select Impact Level",
        choices = choices,
        selected = selected_choice,
        multiple = TRUE,
        width = "100%",
        options = list(placeholder = placeholder_text)
      )
    })

    genotype_results_ui <- function(ns) {
      bslib::card(
        div(
          style = "overflow-x: auto;",
          reactable::reactableOutput(ns("genotype_results_tbl"))
        ),
        class = "mt-3",
        bslib::card_footer(
          textInput(
            ns("File_name"),
            "Enter File Name",
            value = "Chrom 05",
            width = "30%"
          ),
          div(
            style = "display: flex; justify-content: space-between; align-items: center; gap: 10px; width: 100%; flex-grow: 1;",
            shinyWidgets::downloadBttn(
              color = "primary",
              style = "unite",
              icon = icon("upload"),
              outputId = ns("download_excel"),
              label = "Export Genotype Matrix as .xlsx"
            )
          )
        )
      )
    }

    observe({
      req(
        rv$connected,
        values$result$chrom,
        values$result$start,
        values$result$end
      )

      c_type <- rv$conn_type
      d_path <- rv$db_path
      chr <- values$result$chrom
      st <- values$result$start
      en <- values$result$end

      p <- future::future(
        {
          if (c_type == "sqlite") {
            temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
            res <- query_db(
              con = temp_con,
              table_name = "variants",
              chrom = chr,
              start = st,
              end = en
            )
            disconnect_local_db(temp_con, quiet = TRUE)
          } else {
            res <- pg_query_db(
              table_name = "variants",
              chrom = chr,
              start = st,
              end = en
            )
          }
          if (nrow(res) > 0) return(res$variant_id) else return(NULL)
        },
        seed = TRUE
      )

      promises::then(
        p,
        onFulfilled = function(var_ids) {
          if (!is.null(var_ids)) {
            updateSelectizeInput(
              session,
              "manual_variant_ids",
              choices = var_ids,
              server = TRUE
            )
          } else {
            updateSelectizeInput(
              session,
              "manual_variant_ids",
              choices = character(0),
              options = list(placeholder = "No variants found in region")
            )
          }
        },
        onRejected = function(err) {
          message("Failed to fetch variant IDs: ", err$message)
        }
      )
    })

    observeEvent(input$get_pcv_btn, {
      values$query_geno_react <- NULL
      values$filtered_pcvs_by_meta <- NULL
      req(rv$connected)

      updateTabsetPanel(session, "param_header", selected = "pcv_tab")
      updateTabsetPanel(
        session,
        inputId = "pcv_nav_id",
        selected = "PCVs for KASP Marker Design"
      )

      # Determine spinner text and check requirements upfront
      if (input$impact_card == "By Impact & AF") {
        req(input$impact_level, values$result, input$af_range)
        spinner_text <- "Filtering by Impact and AF..."
      } else if (input$impact_card == "By Variant ID") {
        req(input$manual_variant_ids)
        spinner_text <- "Extracting Selected Variants..."
      } else {
        return() # Should not happen
      }

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = spinner_text
      )

      # Capture reactive inputs for the future
      c_type <- rv$conn_type
      d_path <- rv$db_path
      impact_card_val <- input$impact_card
      impact_lvl <- input$impact_level
      res <- values$result
      af_val <- input$af_range
      manual_ids <- input$manual_variant_ids
      meta_cols <- c("chrom", "pos", "ref", "alt", "variant_type", 'major_allele', "minor_allele", "major_allele_freq", "minor_allele_freq")

      p <- future::future({
        # This code runs in a separate process
        if (c_type == "sqlite") {
          temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
          on.exit(disconnect_local_db(temp_con, quiet = TRUE))
        }

        if (impact_card_val == "By Impact & AF") {
          impact_result <- if (c_type == "sqlite") {
            query_by_impact(con = temp_con, impact_level = impact_lvl, chrom = res$chrom, start = res$start, end = res$end)
          } else {
            pg_query_by_impact(impact_level = impact_lvl, chrom = res$chrom, start = res$start, end = res$end)
          }
          if (is.null(impact_result) || nrow(impact_result) == 0) return(list(data = NULL, type = "impact"))

          full_gt <- if (c_type == "sqlite") {
            query_genotypes(con = temp_con, variant_ids = impact_result$variant_id, meta_data = meta_cols)
          } else {
            pg_query_genotypes(variant_ids = impact_result$variant_id, meta_data = meta_cols)
          }
          if (is.null(full_gt) || nrow(full_gt) == 0) return(list(data = NULL, type = "impact"))

          filtered_af_result <- if (c_type == "sqlite") {
            filter_by_af(gt = full_gt, min_af = af_val)
          } else {
            pg_filter_by_af(gt = full_gt, min_af = af_val)
          }
          if (is.null(filtered_af_result) || nrow(filtered_af_result) == 0) return(list(data = NULL, type = "impact"))

          filtered_ids <- filtered_af_result$variant_id
          query_result <- full_gt[full_gt$variant_id %in% filtered_ids, ]
          return(list(data = query_result, type = "impact"))

        } else if (impact_card_val == "By Variant ID") {
          query_result <- if (c_type == "sqlite") {
            query_genotypes(con = temp_con, variant_ids = manual_ids, meta_data = meta_cols)
          } else {
            pg_query_genotypes(variant_ids = manual_ids, meta_data = meta_cols)
          }
          return(list(data = query_result, type = "id"))
        }
      }, seed = TRUE, packages = c("panGenomeBreedr"))

      promises::then(
        p,
        onFulfilled = function(result) {
          values$query_geno_react <- result$data
          if (!is.null(result$data) && nrow(result$data) > 0) {
            if (result$type == "impact") {
              show_toast_success(text = paste("Found", nrow(result$data), "Putative Causal Variants"))
            } else {
              show_toast_success(text = paste("Extracted", nrow(result$data), "Selected Variants"))
            }
          } else {
            if (result$type == "impact") {
              shinyWidgets::show_alert(
                title = "No Variants Found",
                text = paste("No variants found with MAF >=", af_val, ". Try a lower threshold."),
                type = "warning", timer = 5000
              )
            } else {
              shinyWidgets::show_alert(
                title = "No Variants Found",
                text = "Could not extract genotypes for the selected IDs.",
                type = "warning"
              )
            }
          }
          shinybusy::remove_modal_spinner()
        },
        onRejected = function(e) {
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_alert(title = "Error", text = e$message, type = "error")
        }
      )

      output$pcvs_kasp_marker_design_result <- renderUI({
        genotype_results_ui(ns)
      })
    })

    # Evaluate target genotype matrix to compile final output matrix
    output$download_excel <- downloadHandler(
      filename = function() {
        paste0(
          "dbquery_",
          gsub(' ', '', input$File_name),
          "_variant_matrix.xlsx"
        )
      },
      content = function(file) {
        active_data <- if (!is.null(values$filtered_pcvs_by_meta)) {
          values$filtered_pcvs_by_meta
        } else {
          values$query_geno_react
        }
        writexl::write_xlsx(active_data, path = file)
      }
    )

    observeEvent(input$go_back, {
      updateTabsetPanel(session, inputId = 'param_header', selected = 'pcv_tab')
      bslib::toggle_sidebar(id = 'db_sidebar', open = 'open')
      update_sidebar_buttons("get_pcv_sidebar_btn")
    })

    # ==========================================================================
    # KASP MARKER DESIGN LOGIC
    # ==========================================================================

    kasp_des.result <- reactiveVal(NULL)
    kasp_des.plot <- reactiveVal(NULL)

    observeEvent(input$modal_run_but, {
      # Isolate target alleles for locus assay development
      active_data <- if (!is.null(values$filtered_pcvs_by_meta)) {
        values$filtered_pcvs_by_meta
      } else {
        values$query_geno_react
      }

      req(
        active_data,
        input$modal_genome_file$datapath,
        input$modal_maf,
        input$modal_marker_ID
      )
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Designing KASP Marker... Please wait."
      )
      list_markers <- list()
      list_plots <- list()
      tryCatch(
        {
          gt_cols <- colnames(active_data)
          id_col <- gt_cols[grep("id", gt_cols, ignore.case = TRUE)[1]]
          chrom_col <- gt_cols[grep("chro", gt_cols, ignore.case = TRUE)[1]]
          pos_col <- gt_cols[grep("pos", gt_cols, ignore.case = TRUE)[1]]
          ref_col <- gt_cols[grep("ref", gt_cols, ignore.case = TRUE)[1]]
          alt_col <- gt_cols[grep("alt", gt_cols, ignore.case = TRUE)[1]]
          geno_start <- 11
          unique_chr <- unique(active_data[[chrom_col]])

          for (marker in input$modal_marker_ID) {
            result_data <- kasp_marker_design(
              vcf_file = NULL,
              gt_df = active_data,
              variant_id_col = id_col,
              chrom_col = chrom_col,
              pos_col = pos_col,
              ref_al_col = ref_col,
              alt_al_col = alt_col,
              geno_start = input$modal_genotype_start_col,
              marker_ID = marker,
              chr = unique_chr,
              genome_file = input$modal_genome_file$datapath,
              plot_file = tempdir(),
              region_name = input$modal_reg_name,
              maf = input$modal_maf,
              save_alignment = FALSE
            )
            list_markers[[marker]] <- result_data$marker_data
            list_plots[[marker]] <- result_data$plot
          }
          kasp_des.result(data.table::rbindlist(list_markers))
          kasp_des.plot(list_plots)

          shinyWidgets::show_alert(
            title = "Success!",
            text = paste(
              length(list_markers),
              "KASP marker(s) and",
              length(list_plots),
              "plot(s) designed successfully"
            ),
            type = "success",
            timer = 5000
          )
        },
        error = function(e) {
          kasp_des.result(NULL)
          shinyWidgets::show_alert(
            title = "Error",
            text = paste("Error: ", e$message),
            type = "error"
          )
        },
        finally = {
          shinybusy::remove_modal_spinner()
        }
      )
    })

    observe({
      req(kasp_des.result(), kasp_des.plot())
      updateSelectizeInput(
        session,
        inputId = "done_markers",
        choices = names(kasp_des.result())
      )
      updateSelectizeInput(
        session,
        inputId = "plot_choice",
        choices = names(kasp_des.plot())
      )
    })

    output$kasp_table <- DT::renderDT({
      req(kasp_des.result())
      DT::datatable(
        data = kasp_des.result(),
        options = list(scrollX = TRUE, pageLength = 10, serverSide = TRUE),
        escape = TRUE
      )
    })

    output$download_table <- downloadHandler(
      filename = function() {
        clean_name <- gsub("[^[:alnum:]_-]", "_", input$file_name)
        paste0(clean_name, input$exten)
      },
      content = function(file) {
        if (input$exten == ".csv") {
          write.csv(kasp_des.result(), file, row.names = FALSE)
        } else if (input$exten == ".xlsx") {
          openxlsx::write.xlsx(kasp_des.result(), file)
        }
      }
    )

    output$plot <- renderPlot({
      req(kasp_des.plot(), input$plot_choice)
      print(kasp_des.plot()[[input$plot_choice]])
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        clean_marker <- gsub("[^[:alnum:]_-]", "_", input$marker_ID)
        paste0("alignment_", clean_marker, ".pdf")
      },
      content = function(file) {
        plots <- kasp_des.plot()
        grDevices::pdf(file, width = 24, height = 9, onefile = TRUE)
        if (is.list(plots)) {
          lapply(plots, print)
        } else {
          print(plots)
        }
        grDevices::dev.off()
      }
    )
  })
}