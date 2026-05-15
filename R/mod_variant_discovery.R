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

    uiOutput(outputId = ns("variant_hotspot_plot")),

    tags$br(),

    # Center the coordinate entry card
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
          uiOutput(ns("pcvs_kasp_marker_design_result"))
        )
      )
    )
  }

  query_database_card <- function(ns) {
    bslib::card(
      class = "shadow p",
      bslib::card_header(
        tags$strong("Extract Database Records"),
        class = 'text-center',
        style = "font-size:18px;"
      ),
      radioButtons(
        width = "100%",
        inputId = ns("query_database"),
        label = "Query Type:",
        choices = c(
          "Query by Coordinates" = "q_entire",
          "Annotation Summary" = "q_annt"
        ),
        selected = character(0)
      ),
      uiOutput(outputId = ns("display_qd_choice")),
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
    bslib::navset_card_underline(
      id = ns('impact_card'),
      title = tags$span(
      class = "fw-bold",
      "Extract Putative Causal Variants: "
    ),
      
      # ── Filtering by Impact levels ──
      bslib::nav_panel(
        "By Impact & AF",
        icon = icon("filter"),
        div(
          class = "p-3",
          selectInput(
            ns("impact_level"),
            "Impact Level",
            choices = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
            width = "100%"
          ),
          numericInput(
            ns("af_range"),
            "Filter by Alternate Allele Frequency",
            min = 0,
            max = 1,
            value = 0.05,
            width = "100%"
          ),
          verbatimTextOutput(ns('alt_freq_range')),
          div(
            style = "display: flex; justify-content: center; margin-top: 15px;",
            actionButton(ns("get_pcv_btn"), tags$b("Extract PCVs"), class = "btn-warning", width = "70%", icon = icon("filter"))
          )
        )
      ),
      
      # ── Manual ID search ──
      bslib::nav_panel(
        "Select by ID",
        icon = icon("list-check"),
        div(
          class = "p-3",
          selectizeInput(
            ns("manual_variant_ids"),
            "Select Variant IDs",
            choices = NULL,
            multiple = TRUE,
            width = "100%",
            options = list(placeholder = "Select variants from region...")
          ),
          div(
            style = "display: flex; justify-content: center; margin-top: 15px;",
            actionButton(ns("get_manual_pcv_btn"), tags$b("Extract Selected"), class = "btn-warning", width = "70%", icon = icon("check"))
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
        tags$h5("Local Database File (.db)", class = "fw-bold"),
        p(
          "Perform offline variant discovery using a local database file.",
          class = "text-muted mb-4"
        ),
        shinyFiles::shinyFilesButton(
          id = ns("connect_btn"),
          label = "Browse for SQLite File",
          icon = icon("folder-open"),
          title = "Select Database File",
          class = "btn-info btn-lg rounded-pill px-4",
          multiple = FALSE,
        )
      ),
      bslib::card(
        class = "text-center shadow-sm p-4",
        style = "border-top: 4px solid #27AE60;",
        tags$div(class = "h1 text-success mb-3", icon("server")),
        tags$h5("Cloud / Remote Database", class = "fw-bold"),
        p(
          "Perform online variant discovery using hosted pangenome resources.",
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
            ),

            actionButton(
              ns("go_back"),
              "Back to Get PCVs",
              icon = icon("arrow-left", class = "me-2"),
              class = "btn-outline-secondary w-100 py-2 fw-bold shadow-sm"
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
#' @importFrom RSQLite dbConnect dbDisconnect
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
      req(values$query_geno_react)
      bslib::toggle_sidebar(id = 'db_sidebar', open = 'closed')
      updateTabsetPanel(inputId = 'param_header', selected = "mark_design")
      updateSelectizeInput(session, inputId = "modal_marker_ID", choices = values$query_geno_react$variant_id, server = TRUE)
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

    #  SQLite Connection
    volumes <- shinyFiles::getVolumes()
    shinyFiles::shinyFileChoose(
      input,
      "connect_btn",
      roots = volumes,
      session = session
    )

    observeEvent(input$connect_btn, {
      shinyjs::hide(id = 'impact_card')
      file_selected <- shinyFiles::parseFilePaths(volumes, input$connect_btn)

      if (nrow(file_selected) > 0) {
        shinybusy::show_modal_spinner(
          spin = "fading-circle",
          color = "#27AE60",
          text = "Connecting to Database... Please wait."
        )
        tryCatch(
          {
            db_path <- as.character(file_selected$datapath)
            if (!file.exists(db_path)) {
              shinyWidgets::show_toast(
                "Error: Database file not found",
                type = "error"
              )
              return()
            }
            if (!is.null(rv$conn) && rv$conn_type == "sqlite") {
              RSQLite::dbDisconnect(rv$conn)
              rv$conn <- NULL
            }

            rv$conn <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
            rv$connected <- TRUE
            rv$conn_type <- "sqlite"
            rv$tables <- list_sqlite_tables(db_path)
            rv$sample_metadata <- get_sample_metadata(db_path = db_path)
            rv$db_path <- db_path

            shinyjs::hide("pre_connection_panel")
            shinyjs::show("active_dashboard_panel")
            shinyWidgets::show_alert(
              title = "Success!",
              text = "SQLite connected successfully",
              type = "success",
              timer = 3000
            )
          },
          error = function(e) {
            shinyWidgets::show_alert(
              title = "Connection Error!",
              text = paste("Failed to connect or read from database:", e$message),
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
            icon("cloud", class = "text-success me-2", style = "font-size: 1.5rem;"),
            tags$b("Connect to Cloud API")
          ),
          size = "m",
          easyClose = TRUE,
          fade = TRUE,
          div(
            class = "mb-4 text-muted",
            "By default, panGB connects to the public ", tags$b("Sorghum Pangenome"), " database. ",
            "If your institution hosts a custom pangenome database for another crop, you can specify the endpoint below."
          ),
          radioButtons(
            ns("api_choice"),
            label = tags$b("Select Database Source:"),
            choices = c(
              "Public Sorghum Pangenome Database" = "default",
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
      removeModal() # Close the modal upon submission
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Connecting to database..."
      )
      tryCatch(
        {
          # Determine the target API URL based on user selection
          target_url <- if (input$api_choice == "custom") {
            trimws(input$api_url)
          } else {
            "http://16.171.142.87:8000"  # Sorghum pangenome API -- will change  from mini
          }
          
          # Set endpoint API
          set_api_url(target_url)
          # Verify API is alive and get table list via wrapper
          rv$tables <- pg_list_tables()

          rv$connected <- TRUE
          rv$sample_metadata <- pg_get_sample_metadata() # Fetch metadata
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

    # 3. Disconnect
    observeEvent(input$disconnect_btn, {
      if (!is.null(rv$conn) && rv$conn_type == "sqlite") {
        RSQLite::dbDisconnect(rv$conn)
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
          RSQLite::dbDisconnect(rv$conn)
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

      #! QC's to ensure that all columnames in metadata is not passed to selectinput
      # Columns to always exclude
      cols_to_always_exclude <- c("lat", "lon", "latitude", "lib", "sample", "plantname", "pinumber", "array_index", "accession_id", "sample_id")

      candidate_cols <- setdiff(colnames(meta), cols_to_always_exclude)

      # Filter columns to find good candidates for categorization
      suitable_cols <- Filter(function(col_name) {
        col_data <- meta[[col_name]]
        valid_data <- na.omit(col_data)
        if (length(valid_data) == 0) return(FALSE)

        num_unique <- length(unique(valid_data))

        # Exclude columns that are likely unique identifiers
        if (num_unique > 0.9 * length(valid_data) && num_unique > 100) {
          return(FALSE)
        }

        # For numeric columns, only include if they have few unique values
        if (is.numeric(col_data)) {
          return(num_unique < 25)
        }

        #  For character/factor columns, include if they don't have too many unique values
        if (is.character(col_data) || is.factor(col_data)) {
          return(num_unique < 100) 
        }

        return(FALSE) 
      }, candidate_cols)

      updateSelectInput(
        session,
        "color_by_col",
        choices = suitable_cols,
        selected = if ("countryorigin" %in% suitable_cols) "countryorigin" else suitable_cols[1]
      )
    })

    # Display map output
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

    # Helper function
    fetch_data <- function(button_id, task_name, data_slot, sqlite_func_name, pg_func_name) {
      observeEvent(input[[button_id]], {
        req(rv$connected)
        
        # Show a non-blocking toast message 
        shinyWidgets::show_toast(
          title = "Fetching Data",
          text = paste("Retrieving", task_name, "in the background..."),
          type = "info",
          timer = 3000,
          position = "bottom-end"
        )
        
        c_type <- rv$conn_type
        d_path <- rv$db_path
        
        # Launch data extraction asynchronously
        p <- future::future({
          if (c_type == "sqlite") {
            do.call(sqlite_func_name, list(db_path = d_path))
          } else {
            do.call(pg_func_name, list())
          }
        }, seed = TRUE, packages = c("panGenomeBreedr"))
        
        # Handle the promise resolution
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
            shinyWidgets::show_alert(title = "Failed!", text = err$message, type = "danger")
          }
        )
      })
    }

    # Fetch Variant Impact Summary
    fetch_data("get_impact", "Variant Impact Summary", "variant_impact",
               "variant_impact_summary", "pg_variant_impact_summary")

    # Fetch Variant Statistics
    fetch_data("get_stats", "Variant Statistics", "variant_stats",
               "variant_stats", "pg_variant_stats")

    # Fetch Table Summaries
    fetch_data("get_summary", "Table Summaries", "sqlite_summary",
               "summarize_sqlite_tables", "pg_summarize_tables")

    # Fetch Variant Type Counts
    fetch_data("get_types", "Variant Type Counts", "variant_count",
               "count_variant_types", "pg_count_variant_types")

    # Fetch Table Schema
    observeEvent(input$get_schema, {
      req(rv$connected, input$table_name_lst)
      
      shinyWidgets::show_toast(
        title = "Fetching Data",
        text = paste("Retrieving schema for", input$table_name_lst, "in the background..."),
        type = "info",
        timer = 3000,
        position = "bottom-end"
      )
      
      c_type <- rv$conn_type
      d_path <- rv$db_path
      t_name <- input$table_name_lst
      
      p <- future::future({
        if (c_type == "sqlite") {
          list_table_columns(db_path = d_path, table_name = t_name)
        } else {
          pg_list_table_columns(table_name = t_name)
        }
      }, seed = TRUE, packages = c("panGenomeBreedr"))
      
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
          shinyWidgets::show_alert(title = "Failed!", text = err$message, type = "danger")
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
      last_action = NULL
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

  

    # Widget for gff3 file setting.
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



    # Widget for manual cordinate setting
    observeEvent(input$set_genocod_btn, {
      req(input$chrom, input$start, input$end)
      tryCatch(
        {
          values$result <- list(
            chrom = sprintf("Chr%02d", input$chrom),
            start = input$start,
            end = input$end
          )
          shinyWidgets::show_alert(
            title = "Gene Coordinates Set",
            text = sprintf(
              "Chromosome: %s | Start: %d | End: %d",
              sprintf("Chr%02d", input$chrom),
              input$start,
              input$end
            ),
            type = "success",
            timer = 5000
          )
          shinyjs::show(id = 'impact_card')
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Error Setting Coordinates",
            text = e$message,
            type = "error"
          )
        }
      )
    })


    # VALUE BOXES UI
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
        
        # Calculate and format the length of the genomic region
        l <- values$result$end - values$result$start + 1
        if (l < 1000) {
          region_len <- paste0(l, " bp")
        } else if (l < 1e6) {
          region_len <- paste0(round(l / 1000, 2), " kb")
        } else {
          region_len <- paste0(round(l / 1e6, 2), " Mb")
        }
      }
      
      # Visual format for value box elements
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


    # Variant Hotspot plot.
    hotspot_result <- reactive({
      req(rv$connected, values$result$chrom, values$result$start, values$result$end)
      
      # Extract reactive variables to local scope to safely pass to the future
      c_type <- rv$conn_type
      d_path <- rv$db_path
      chr <- values$result$chrom
      st <- values$result$start
      en <- values$result$end
      
      # Launch data extraction and plotting asynchronously
      future::future({
        if (c_type == "sqlite") {
          v_table <- query_db(
            db_path = d_path,
            table_name = "variants",
            chrom = chr,
            start = st,
            end = en
          )
          a_table <- query_db(
            db_path = d_path,
            table_name = "annotations",
            chrom = chr,
            start = st,
            end = en
          )
        } else {
          v_table <- pg_query_db(
            table_name = "variants",
            chrom = chr,
            start = st,
            end = en
          )
          a_table <- pg_query_db(
            table_name = "annotations",
            chrom = chr,
            start = st,
            end = en
          )
        }
        
        plot_variant_hotspots(
          variant_table = v_table,
          annotation_table = a_table,
          region_start = st,
          region_end = en
        )
      }, seed = TRUE, packages = c("ggplot2", "scales"))
    })
    
    output$variant_hotspot_plot <- renderUI({
      req(values$result) # Only depend on coordinates being set so the UI renders instantly
      bslib::card(
        class = "shadow-sm border-0",
        bslib::card_header(
          class = "bg-light d-flex align-items-center",
          icon("chart-line", class = "me-2 text-primary"),
          tags$strong("Variant Hotspot Plot")
        ),
        bslib::card_body(
          # Show a spinner while the promise from the future resolves
          shinycssloaders::withSpinner(
            plotOutput(outputId = ns('hotspot_plot')),
            type = 4, color = "#27AE60",
            caption = "Plotting Variant Hotspots..."
          )
        )
      )
    })

    output$hotspot_plot <- renderPlot({
      req(hotspot_result())
      # Return the promise directly
      hotspot_result()
    })
    
    

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
      updateTabsetPanel(
        session,
        inputId = "query_db_nav_id",
        selected = "Main Database Results"
      )
      values$last_action <- NULL
      req(rv$connected, values$result, input$query_database)

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Querying Database... Please wait."
      )

      tryCatch(
        {
          if (!is.null(input$table_name) && input$table_name != "") {
            # Branching Based on DB Type
            if (rv$conn_type == "sqlite") {
              values$query_db_val <- query_db(
                db_path = rv$db_path,
                table_name = input$table_name,
                chrom = values$result$chrom,
                start = values$result$start,
                end = values$result$end,
                gene_name = if (input$query_gene_name == "") {
                  NULL
                } else {
                  input$query_gene_name
                }
              )
            } else if (rv$conn_type == "postgres") {
              values$query_db_val <- pg_query_db(
                table_name = input$table_name,
                chrom = values$result$chrom,
                start = values$result$start,
                end = values$result$end,
                gene_name = if (input$query_gene_name == "") {
                  NULL
                } else {
                  input$query_gene_name
                }
              )
            }
            values$last_action <- "main_db"
          }

          if (!is.null(input$table_name_v) && !is.null(input$table_name_a)) {
            # Branching Based on DB Type
            if (rv$conn_type == "sqlite") {
              values$query_ann_react <- query_ann_summary(
                db_path = rv$db_path,
                variants_table = input$table_name_v,
                annotations_table = input$table_name_a,
                chrom = values$result$chrom,
                start = values$result$start,
                end = values$result$end
              )
            } else if (rv$conn_type == "postgres") {
              # Assumes API takes the same params. Adjust if pg_ requires different arguments.
              values$query_ann_react <- pg_query_ann_summary(
                annotations_table = input$table_name_a,
                variants_table = input$table_name_v,
                chrom = values$result$chrom,
                start = values$result$start,
                end = values$result$end
              )
            }
            values$last_action <- "annotation_summ"
          }
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
          shinyjs::delay(ms = 1000, {
            shinybusy::remove_modal_spinner()
            if (!is.null(values$last_action)) {
              if (values$last_action == "annotation_summ") {
                show_toast_success("Variant Filtered By Annotation Summary")
              } else if (values$last_action == "main_db") {
                show_toast_success(text = paste("Queried by", input$table_name))
              }
            }
          })
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

    observe({
      req(input$query_database)
      if (
        !is.null(input$query_database) && input$query_database == "q_entire"
      ) {
        output$query_db_display <- renderUI({
          reactable::renderReactable({
            req(values$query_db_val)
            render_reactable(values$query_db_val)
          })
        })
      } else if (
        !is.null(input$query_database) && input$query_database == "q_annt"
      ) {
        output$query_db_display <- renderUI({
          annotation_summary_results_ui(ns)
        })
      }
    })

    # New observer for the sidebar button to get PCVs
    observeEvent(input$get_pcv_sidebar_btn, {
      updateTabsetPanel(session, "param_header", selected = "pcv_tab")
      update_sidebar_buttons("get_pcv_sidebar_btn")
    })


    # ==========================================================================
    # CAUSAL VARIANTS (PCVs) EXTRACTION
    # ==========================================================================
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
            # ,
            # shinyWidgets::actionBttn(
            #   style = "unite",
            #   color = "success",
            #   inputId = ns("push_1"),
            #   label = "Design KASP Markers for PCVs",
            #   icon = icon("dna")
            # )
          )
        )
      )
    }

       # Background fetching of available variants for the manual dropdown
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
            res <- query_db(
              db_path = d_path,
              table_name = "variants",
              chrom = chr,
              start = st,
              end = en
            )
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

    # Execute query for manually selected variant IDs
    observeEvent(input$get_manual_pcv_btn, {
      req(input$manual_variant_ids, rv$connected)

      values$query_geno_react <- NULL
      updateTabsetPanel(session, "param_header", selected = "pcv_tab")
      updateTabsetPanel(
        session,
        "pcv_nav_id",
        selected = "PCVs for KASP Marker Design"
      )

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Extracting Selected Variants... Please wait."
      )

      tryCatch(
        {
          if (rv$conn_type == "sqlite") {
            values$query_geno_react <- query_genotypes(
              db_path = rv$db_path,
              variant_ids = input$manual_variant_ids,
              meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
            )
          } else {
            values$query_geno_react <- pg_query_genotypes(
              variant_ids = input$manual_variant_ids,
              meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
            )
          }

          if (
            !is.null(values$query_geno_react) &&
              nrow(values$query_geno_react) > 0
          ) {
            show_toast_success(
              text = paste(
                "Extracted",
                nrow(values$query_geno_react),
                "Selected Variants"
              )
            )
          } else {
            shinyWidgets::show_alert(
              title = "No Variants Found",
              text = "Could not extract genotypes for the selected IDs.",
              type = "warning"
            )
          }
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Error",
            text = e$message,
            type = "error"
          )
        },
        finally = {
          shinybusy::remove_modal_spinner()
        }
      )

      output$pcvs_kasp_marker_design_result <- renderUI({
        genotype_results_ui(ns)
      })
    })

    # Branching Based on DB Type
    query_by_impact_result <- reactive({
      req(rv$connected, input$impact_level, values$result)
      if (rv$conn_type == "sqlite") {
        query_by_impact(
          db_path = rv$db_path,
          impact_level = input$impact_level,
          chrom = values$result$chrom,
          start = values$result$start,
          end = values$result$end
        )
      } else {
         pg_query_by_impact(
          impact_level = input$impact_level,
          chrom = values$result$chrom,
          start = values$result$start,
          end = values$result$end
        )
      }
    })

    # Branching Based on DB Type
    hold_genotypes_impact <- reactive({
      req(query_by_impact_result(), rv$connected)
      if (rv$conn_type == "sqlite") {
        query_genotypes(
          db_path = rv$db_path,
          variant_ids = query_by_impact_result()$variant_id,
          meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
        )
      } else {
         pg_query_genotypes(
          variant_ids = query_by_impact_result()$variant_id,
          meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
        )
      }
    })


    calc_af_result <- reactive({
      # Require the genotype matrix (cached from either SQLite or PG)
      gt_data <- hold_genotypes_impact()
      req(gt_data)
      # Calc aternate allele frequency
      tryCatch({
        if (rv$conn_type == "sqlite") {
          alt_af_df <- calc_af(
            gt = gt_data,
            variant_id_col = 'variant_id',
            chrom_col = 'chrom',
            pos_col = 'pos'
          )
        } else {
          alt_af_df <- pg_calc_af(
            gt = gt_data,
            variant_id_col = 'variant_id',
            chrom_col = 'chrom',
            pos_col = 'pos'
          )
        }
        
        # Extract the range for the UI summary
        alt_af_range <- range(alt_af_df$alt_af, na.rm = TRUE)

        if (all(!is.finite(alt_af_range))) {
          return(NULL)
        }
        
        return(alt_af_range)
        
      }, error = function(e) {
        message("Error in AF calculation: ", e$message)
        return(NULL)
      })
    })


    output$alt_freq_range <- renderPrint({
      req(calc_af_result())
      alt_af_range <- calc_af_result()
      if (is.null(alt_af_range) || length(alt_af_range) != 2) {
        cat(
          "Note: Alternate allele frequency for this impact level is unavailable."
        )
      } else {
        cat(
          "Alternate allele frequency range:",
          alt_af_range[1],
          "-",
          alt_af_range[2]
        )
      }
    })

    # Filter Logic - Assuming local computation is fine once fetched, or pg_filter
    query_by_alf_result <- reactive({
      req(input$af_range, hold_genotypes_impact())
      if (rv$conn_type == "sqlite") {
        filter_by_af(gt = hold_genotypes_impact(), min_af = input$af_range)
      } else {
        # Assuming you pipe the genotypes into the pg wrapper as done in your test script
         pg_filter_by_af(
          gt = hold_genotypes_impact(),
          min_af = input$af_range
        )
      }
    })

    observeEvent(input$get_pcv_btn, {
      values$query_geno_react <- NULL
      req(query_by_alf_result(), rv$connected)
      
      updateTabsetPanel(session, "param_header", selected = "pcv_tab")
      
      updateTabsetPanel(
        session,
        inputId = "pcv_nav_id",
        selected = "PCVs for KASP Marker Design"
      )

      tryCatch(
        {
          # Eliminate the redundant database query by simply subsetting the data already in memory
          full_gt <- hold_genotypes_impact()
          filtered_ids <- query_by_alf_result()$variant_id
          values$query_geno_react <- full_gt[full_gt$variant_id %in% filtered_ids, ]
          
          if (!is.null(values$query_geno_react) && nrow(values$query_geno_react) > 0) {
            show_toast_success(
              text = paste(
                "Found",
                nrow(values$query_geno_react),
                "Putative Causal Variants"
              )
            )
          } else {
            shinyWidgets::show_alert(
              title = "No Variants Found",
              text = paste(
                "No putative causal variants found at MAF threshold of",
                input$af_range,
                ". Try lowering the MAF threshold."
              ),
              type = "warning",
              timer = 5000
            )
          }
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Error",
            text = e$message,
            type = "error"
          )
        }
      )
      output$pcvs_kasp_marker_design_result <- renderUI({
        genotype_results_ui(ns)
      })
    })

    output$genotype_results_tbl <- reactable::renderReactable({
      req(values$query_geno_react)
      render_reactable(values$query_geno_react)
    })

    output$download_excel <- downloadHandler(
      filename = function() {
        paste0(
          "dbquery_",
          gsub(' ', '', input$File_name),
          "_variant_matrix.xlsx"
        )
      },
      content = function(file) {
        writexl::write_xlsx(values$query_geno_react, path = file)
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
      req(
        values$query_geno_react,
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
          gt_cols <- colnames(values$query_geno_react)
          id_col <- gt_cols[grep("id", gt_cols, ignore.case = TRUE)[1]]
          chrom_col <- gt_cols[grep("chro", gt_cols, ignore.case = TRUE)[1]]
          pos_col <- gt_cols[grep("pos", gt_cols, ignore.case = TRUE)[1]]
          ref_col <- gt_cols[grep("ref", gt_cols, ignore.case = TRUE)[1]]
          alt_col <- gt_cols[grep("alt", gt_cols, ignore.case = TRUE)[1]]
          geno_start <- 7
          unique_chr <- unique(values$query_geno_react[[chrom_col]])

          for (marker in input$modal_marker_ID) {
            result_data <- kasp_marker_design(
              vcf_file = NULL,
              gt_df = values$query_geno_react,
              variant_id_col = id_col,
              chrom_col = chrom_col,
              pos_col = pos_col,
              ref_al_col = ref_col,
              alt_al_col = alt_col,
              geno_start = geno_start,
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