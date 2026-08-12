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
          "Geographic Distribution of Accessions"
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
              "Database Summary"
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
              "Table Schema"
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
              "Select Table",
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
  bslib::card(
    class = "shadow-sm",
    bslib::card_header(
      tags$h5("Step 1: Define Genomic Region", class = "fw-bold m-0")
    ),
    bslib::card_body(
      radioButtons(
        ns("define_region_method"),
        label = "Method:",
        choices = c("By Gene ID" = "gene", "By Coordinates" = "coords"),
        selected = "gene",
        inline = TRUE
      ),
      tags$hr(),
      # Conditional Panel for Gene ID
      conditionalPanel(
        condition = paste0("input['", ns("define_region_method"), "'] === 'gene'"),
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
        radioButtons(
          ns("input_method"),
          tags$span(
            class = "text-muted small text-uppercase fw-bold",
            "GFF3 File Source"
          ),
          choices = c("Remote URL" = "url", "Local Upload" = "file"),
          selected = "url",
          inline = TRUE
        ),
        conditionalPanel(
          condition = paste0("input['", ns("input_method"), "'] === 'url'"),
          textInput(
            ns("gff_url"),
            NULL,
            width = "100%",
            value = "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz"
          )
        ),
        conditionalPanel(
          condition = paste0("input['", ns("input_method"), "'] === 'file'"),
          fileInput(
            ns("gff_file"),
            NULL,
            accept = c(".gff3", ".gff", ".gff3.gz", ".gff.gz", "application/gzip", "application/x-gzip"),
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
      ),
      # Conditional Panel for Coordinates
      conditionalPanel(
        condition = paste0("input['", ns("define_region_method"), "'] === 'coords'"),
        numericInput(
          ns("chrom"),
          tags$span(
            class = "text-muted small text-uppercase fw-bold",
            "Chromosome"
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
            "Set Region",
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
    # Step 1: Define Region
    bslib::layout_columns(
      col_widths = c(-2, 8, -2),
      coordinate_entry_card(ns)
    ),

    # Display the selected genomic range
    uiOutput(ns("genomic_range_vboxes")),

    # Step 2: Explore Variants (UI is rendered from server)
    uiOutput(ns("step2_explore_ui"))
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
        "Identify Variants",
        icon = icon("crosshairs"),
        width = "100%",
        class = "btn-outline-primary btn-lg mb-4 fw-bold text-start"
      ),

      # Design KASP Markers button
      actionButton(
        ns("design_kasp_sidebar_btn"),
        "Design Markers",
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
          tags$h4("KASP Design KASP Markers", class = "m-0 fw-bold text-primary")
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
                  "First Genotype Column",
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
                        "Download Format",
                        choices = c(".csv", ".xlsx"),
                        selected = ".xlsx"
                      ),
                      textInput(
                        ns("file_name"),
                        "File Prefix",
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

    # Helper function to apply a list of filters to a dataframe
    apply_filters_to_df <- function(df, filters) {
      if (is.null(df) || is.null(filters) || length(filters) == 0) {
        return(df)
      }
      filtered_df <- df
      for (f in filters) {
        if (!f$column %in% names(filtered_df)) next

        # Drop rows with NA in the filter column to avoid errors
        filtered_df <- filtered_df[!is.na(filtered_df[[f$column]]), ]

        if (is.numeric(f$value) && is.numeric(filtered_df[[f$column]])) {
          filtered_df <- filtered_df[filtered_df[[f$column]] >= f$value[1] & filtered_df[[f$column]] <= f$value[2], ]
        } else {
          filtered_df <- filtered_df[as.character(filtered_df[[f$column]]) %in% as.character(f$value), ]
        }
      }
      return(filtered_df)
    }

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
            rv$tables <- list_tables(con = rv$conn)
            rv$sample_metadata <- fetch_accession_metadata(con = rv$conn)
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
          rv$tables <- list_tables(connect_db_mode = 'online')
          rv$connected <- TRUE
          rv$sample_metadata <- fetch_accession_metadata(connect_db_mode = 'online')
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
      req(rv$sample_metadata, input$color_by_col)
      plot_accession_map(
        metadata = rv$sample_metadata, color_by = input$color_by_col
      )
    })

    # ------------------------------------------------------------------
    # DATABASE INFO TAB
    # ------------------------------------------------------------------

    # Refactored fetch_data to work with unified functions that handle both
    # local and online modes, reducing code duplication.
    fetch_data <- function(button_id, task_name, data_slot, unified_func, ...) {
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
        # Capture additional arguments for the function call
        extra_args <- list(...)

        p <- future::future({
          # Prepare arguments based on connection type
          args <- if (c_type == "sqlite") {
            temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
            on.exit(disconnect_local_db(temp_con, quiet = TRUE))
            c(list(con = temp_con), extra_args)
          } else { # postgres
            c(list(connect_db_mode = 'online'), extra_args)
          }
          
          # Call the unified function with the correct arguments
          do.call(unified_func, args)
        }, seed = TRUE, packages = c("panGenomeBreedr"))

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
      summarize_variant_impacts
    )

    fetch_data(
      "get_stats",
      "Variant Statistics",
      "variant_stats",
      summarize_variants,
      include_annotations = FALSE # Explicitly set to FALSE as in original logic
    )

    fetch_data(
      "get_summary",
      "Table Summaries",
      "sqlite_summary",
      summarize_database
    )

    fetch_data(
      "get_types",
      "Variant Type Counts",
      "variant_count",
      count_variant_types
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
            res <- list_columns(con = temp_con, table_name = t_name)
            disconnect_local_db(temp_con, quiet = TRUE)
            res
          } else {
            list_columns(table_name = t_name, connect_db_mode = 'online')
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
      last_action = NULL,
      data_for_marker_design = NULL,
      # New reactive values for on-demand tables
      table_genotypes = NULL,
      table_variants = NULL,
      table_annotations = NULL,
      table_annotation_summary = NULL,
      ld_results = NULL,
      filters_genotypes = NULL,
      filters_variants = NULL,
      filters_annotations = NULL,
      filtered_ids_genotypes = NULL,
      filtered_ids_variants = NULL,
      filtered_ids_annotations = NULL,
      final_selected_variants = NULL,
      plot_filters_annotations = NULL, # For de-coupling plot reactivity
      plot_filters_genotypes = NULL,
      plot_filters_variants = NULL,
      filtered_sample_ids = NULL,
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

    observeEvent(values$result, {
      req(rv$connected, values$result)

      values$ld_results <- NULL
      values$geodesic_plot_obj <- NULL
      # Clear on-demand tables when region changes
      values$table_genotypes <- NULL
      values$table_variants <- NULL
      values$table_annotations <- NULL
      values$table_annotation_summary <- NULL
      values$final_selected_variants <- NULL
      values$filters_genotypes <- NULL; values$filters_variants <- NULL; values$filters_annotations <- NULL
      values$plot_filters_annotations <- NULL; values$plot_filters_genotypes <- NULL; values$plot_filters_variants <- NULL
      filter_states$genotypes <- NULL; filter_states$variants <- NULL; filter_states$annotations <- NULL
      values$filtered_sample_ids <- NULL

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Loading regional data..."
      )

      c_type <- rv$conn_type
      d_path <- rv$db_path
      chr <- values$result$chrom
      st <- values$result$start
      en <- values$result$end

      p <- future::future(
        {
          if (c_type == "sqlite") {
            temp_con <- connect_local_db(d_path, quiet = TRUE)
            on.exit(disconnect_local_db(temp_con, quiet = TRUE))
            list(
              variants = fetch_table_region(con = temp_con, table_name = "variants", chrom = chr, start = st, end = en),
              annotations = fetch_table_region(con = temp_con, table_name = "annotations", chrom = chr, start = st, end = en),
              genotypes = fetch_table_region(con = temp_con, table_name = "genotypes", chrom = chr, start = st, end = en),
              summary = summarize_annotations(con = temp_con, chrom = chr, start = st, end = en)
            )
          } else {
            list(
              variants = fetch_table_region(connect_db_mode = 'online', table_name = "variants", chrom = chr, start = st, end = en),
              annotations = fetch_table_region(connect_db_mode = 'online', table_name = "annotations", chrom = chr, start = st, end = en),
              genotypes = fetch_table_region(connect_db_mode = 'online', table_name = "genotypes", chrom = chr, start = st, end = en),
              summary = summarize_annotations(connect_db_mode = 'online', chrom = chr, start = st, end = en)
            )
          }
        },
        seed = TRUE,
        packages = c("panGenomeBreedr")
      )

      promises::then(
        p,
        onFulfilled = function(data) {
          values$table_variants <- data$variants
          values$table_annotations <- data$annotations
          values$table_genotypes <- data$genotypes
          values$table_annotation_summary <- data$summary
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_toast("Regional data loaded and plot rendered.", type = "success")
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
          gff_df <- gene_coord_gff(trimws(input$gene_name), gff_path)

          # Extract coordinates for the 'gene' feature to standardize values$result
          gene_feature <- gff_df[gff_df$Feature == "gene", ]

          if (nrow(gene_feature) == 0) {
            stop("Could not find 'gene' feature in the GFF file for the provided gene name.")
          }

          # Standardize values$result to be a list with the overall gene boundaries
          values$result <- list(
            chrom = gene_feature$Chromosome[1],
            start = min(gene_feature$Start, na.rm = TRUE),
            end = max(gene_feature$End, na.rm = TRUE)
          )

          shinyWidgets::show_alert(
            title = "Found Gene Co-ordinates",
            text = sprintf("Chromosome: %s | Start: %s | End: %s",
                           values$result$chrom,
                           format(values$result$start, big.mark = ","),
                           format(values$result$end, big.mark = ",")),
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
              res <- summarize_variants(con = temp_con, include_annotations = FALSE)
              disconnect_local_db(temp_con, quiet = TRUE)
              res
            } else {
              summarize_variants(connect_db_mode = 'online', include_annotations = FALSE)
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
              text = sprintf("Chromosome: %s | Start: %s | End: %s",
                             chr_val, format(st_val, big.mark = ","),
                             format(en_val, big.mark = ",")),
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

    output$step2_explore_ui <- renderUI({
      req(values$result)
      tagList(
        tags$hr(class = "my-4"),
        tags$h5("Step 2: Explore Regional Variants", class = "fw-bold mb-3"),

        # Central Filtering Controls
        bslib::card(
          class = "shadow-sm mb-4",
          bslib::card_body(
            class = "p-2",
            div(
              class = "d-flex justify-content-between align-items-center",
              tags$strong(textOutput(ns("filtered_variant_summary_text"))),
              div(
                class = "d-flex gap-2",
                actionButton(ns("refine_plot_btn"), "Apply Filters & Update Plot", icon = icon("search-plus"), class = "btn-primary"),
                actionButton(ns("design_markers_from_region_btn"), "Design Markers", icon = icon("dna"), class = "btn-success"),
                actionButton(ns("clear_all_filters_btn"), "Clear All Filters", icon = icon("xmark"), class = "btn-outline-secondary")
              )
            )
          )
        ),

        # Variant Hotspot Plot is the main feature
        uiOutput(outputId = ns("variant_hotspot_plot_ui")),

        # Accordion for detailed data tables
        bslib::accordion(
          id = ns("explore_data_accordion"),
          open = FALSE,
          multiple = TRUE,
          class = "mt-3",

          # Genotypes Panel
          bslib::accordion_panel(
            title = "Genotypes",
            icon = icon("table"),
            div(
              class = "d-flex justify-content-between align-items-center mb-3",
              tags$strong(textOutput(ns("filtered_sample_count"))),
              actionButton(ns("filter_genotypes_btn"), "Filter Table", icon = icon("filter"), class = "btn-outline-secondary")
            ),
            uiOutput(ns("regional_genotypes_ui"))
          ),

          # Variants Panel
          bslib::accordion_panel(
            title = "Variants",
            icon = icon("table"),
            div(
              class = "d-grid gap-2 d-md-flex justify-content-md-end",
              actionButton(ns("filter_variants_btn"), "Filter Table", icon = icon("filter"), class = "btn-outline-secondary mb-3")
            ),
            uiOutput(ns("regional_variants_ui"))
          ),

          # Annotations Panel
          bslib::accordion_panel(
            title = "Annotations",
            icon = icon("table"),
            div(
              class = "d-grid gap-2 d-md-flex justify-content-md-end",
              actionButton(ns("filter_annotations_btn"), "Filter Table", icon = icon("filter"), class = "btn-outline-secondary mb-3")
            ),
            uiOutput(ns("regional_annotations_ui")),
            tags$hr(),
            bslib::accordion(
              id = ns("annotation_summary_accordion"),
              open = FALSE,
              bslib::accordion_panel(
                "Annotation Summary",
                icon = icon("bars-staggered"),
                uiOutput(ns("annotation_summary_ui"))
              )
            )
          ),

          # LD Analysis Panel
          bslib::accordion_panel(
            title = "Linkage Disequilibrium Analysis",
            icon = icon("link"),
            uiOutput(ns("ld_analysis_ui"))
          )
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
      req(values$result, values$table_annotations, values$table_genotypes)

      # Initialize with full regional data
      plot_ann <- values$table_annotations
      plot_geno <- values$table_genotypes

      # Apply the final intersection of variants from the main filter button.
      if (!is.null(values$final_selected_variants)) {
        plot_ann <- plot_ann[plot_ann$variant_id %in% values$final_selected_variants, ]
        plot_geno <- plot_geno[plot_geno$variant_id %in% values$final_selected_variants, ]
      }

      # Re-apply annotation filters to prevent plotting unfiltered impacts for selected variants.
      if (!is.null(values$plot_filters_annotations)) {
        plot_ann <- apply_filters_to_df(plot_ann, values$plot_filters_annotations)
      }

      # Deduplicate annotations, keeping only the most severe impact per variant.
      if (nrow(plot_ann) > 0 && "impact" %in% names(plot_ann)) {
        impact_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN")
        plot_ann$impact_rank <- factor(plot_ann$impact, levels = impact_levels, ordered = TRUE)
        plot_ann <- plot_ann[order(plot_ann$impact_rank, na.last = TRUE), ]
        plot_ann <- plot_ann[!duplicated(plot_ann$variant_id), ]
        plot_ann$impact_rank <- NULL
      }

      # Handle empty plot case.
      if (nrow(plot_ann) == 0 || nrow(plot_geno) == 0) {
        empty_plot <- ggplot2::ggplot() + ggplot2::theme_void() +
          ggplot2::labs(title = "No variants to display based on current filters.") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        if (input$define_region_method == 'gene') {
          req(input$gene_name)
          gff_path <- if (input$input_method == "url") req(input$gff_url) else req(input$gff_file$datapath)
          gene_model_plot <- plot_gene_model(gene_df = gene_coord_gff(gene_name = input$gene_name, gff_path))
          return(empty_plot / gene_model_plot + patchwork::plot_layout(heights = c(4, 1)))
        } else {
          return(empty_plot)
        }
      }

      # Render plots with the clean, deduplicated data.
      if (input$define_region_method == 'gene') {
        req(input$gene_name)
        gff_path <- if (input$input_method == "url") req(input$gff_url) else req(input$gff_file$datapath)
        hotspot_overlay_plot(
          gene_name = input$gene_name,
          gff_path = gff_path,
          annotations_df = plot_ann,
          genotypes_df = plot_geno,
          selected_variants = unique(values$final_selected_variants)
        )
      } else { # 'coords'
        var_df <- merge(plot_ann, plot_geno, by = "variant_id")
        plot_variant_hotspot(var_df = var_df, start_pos = values$result$start, end_pos = values$result$end)
      }
    })

    output$hotspot_plot <- renderPlot({
      hotspot_plot_obj()
    })

    # ==========================================================================
    # INTERACTIVE FILTERING LOGIC
    # ==========================================================================

    filter_states <- reactiveValues(genotypes = NULL, variants = NULL, annotations = NULL, sample_filters = NULL)
    
    show_filter_modal <- function(table_name, data) {
      req(data)
      # filterable_cols <- names(data)[sapply(data, function(x) !is.list(x) && !is.data.frame(x))]

      showModal(modalDialog(
        title = paste("Filter", tools::toTitleCase(table_name), "Table"),
        size = "l", easyClose = TRUE,
        div(
          class = "p-3 mb-3 bg-light rounded border",
          h5("Add New Filter"),
          bslib::layout_columns(
            col_widths = c(4, 5, 3),
            selectizeInput(ns(paste0("filter_col_", table_name)), "Column", choices = NULL),
            uiOutput(ns(paste0("filter_widget_", table_name))),
            div(class = "d-flex align-items-end", actionButton(ns(paste0("add_filter_", table_name)), "Add Filter", icon = icon("plus"), class = "btn-primary w-100"))
          )
        ),
        h5("Active Filters"),
        uiOutput(ns(paste0("active_filters_ui_", table_name))),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns(paste0("clear_modal_filters_", table_name)), "Clear Filters in Modal", class = "btn-outline-danger"),
          actionButton(ns(paste0("apply_modal_filters_", table_name)), "Apply & Close", icon = icon("check"), class = "btn-success")
        )
      ))
    }
    
    handle_table_filtering <- function(table_name, data_reactive) {
      observeEvent(input[[paste0("filter_", table_name, "_btn")]], {
        req(data_reactive())
        df <- data_reactive()
        show_filter_modal(table_name, df)

        # Use server-side selectize for performance.
        filterable_cols <- names(df)[sapply(df, function(x) !is.list(x) && !is.data.frame(x))]
        updateSelectizeInput(
          session,
          paste0("filter_col_", table_name),
          choices = filterable_cols,
          server = TRUE
        )
      })
    
      output[[paste0("filter_widget_", table_name)]] <- renderUI({
        req(input[[paste0("filter_col_", table_name)]])
        col <- input[[paste0("filter_col_", table_name)]]
        col_data <- data_reactive()[[col]]
        if (is.numeric(col_data)) {
          rng <- range(col_data, na.rm = TRUE)
          sliderInput(ns(paste0("filter_val_", table_name)), "Range", min = rng[1], max = rng[2], value = rng)
        } else {
          choices <- sort(unique(na.omit(as.character(col_data))))
          selectizeInput(ns(paste0("filter_val_", table_name)), "Values", choices = choices, multiple = TRUE, options = list(placeholder = "Select value(s)..."))
        }
      })
    
      observeEvent(input[[paste0("add_filter_", table_name)]], {
        req(input[[paste0("filter_col_", table_name)]], input[[paste0("filter_val_", table_name)]])
        new_filter <- list(column = input[[paste0("filter_col_", table_name)]], value = input[[paste0("filter_val_", table_name)]])
        current_filters <- filter_states[[table_name]]
        filter_states[[table_name]] <- if (is.null(current_filters)) list(new_filter) else c(current_filters, list(new_filter))
      })
    
      output[[paste0("active_filters_ui_", table_name)]] <- renderUI({
        filters <- filter_states[[table_name]]
        if (is.null(filters) || length(filters) == 0) return(tags$p("No filters added yet.", class = "text-muted"))
        tagList(lapply(seq_along(filters), function(i) {
          f <- filters[[i]]; val_str <- paste(f$value, collapse = ", "); if (is.numeric(f$value)) val_str <- paste(f$value[1], "-", f$value[2])
          div(class = "d-flex justify-content-between align-items-center p-2 border-bottom",
              tags$span(tags$b(f$column), ": ", val_str),
              actionButton(ns(paste0("remove_filter_", table_name, "_", i)), "", icon = icon("xmark"), class = "btn-sm btn-outline-danger"))
        }))
      })
    
      observe({
        filters <- filter_states[[table_name]]
        if (!is.null(filters) && length(filters) > 0) {
          lapply(seq_along(filters), function(i) {
            observeEvent(input[[paste0("remove_filter_", table_name, "_", i)]], {
              current_filters <- filter_states[[table_name]]; current_filters[[i]] <- NULL
              filter_states[[table_name]] <- Filter(Negate(is.null), current_filters)
            }, ignoreInit = TRUE, once = TRUE)
          })
        }
      })
    
      observeEvent(input[[paste0("clear_modal_filters_", table_name)]], { filter_states[[table_name]] <- NULL })
    
      observeEvent(input[[paste0("apply_modal_filters_", table_name)]], {
        current_col <- input[[paste0("filter_col_", table_name)]]
        current_val <- input[[paste0("filter_val_", table_name)]]
    
        if (!is.null(current_col) && nzchar(current_col) && !is.null(current_val) && length(current_val) > 0) {
          new_filter <- list(column = current_col, value = current_val)
          is_duplicate <- any(sapply(filter_states[[table_name]], function(f) isTRUE(all.equal(f, new_filter))))
          if (!is_duplicate) {
            filter_states[[table_name]] <- c(filter_states[[table_name]], list(new_filter))
          }
        }
    
        modal_filters <- filter_states[[table_name]]
        values[[paste0("filters_", table_name)]] <- modal_filters
    
        if (!is.null(modal_filters) && length(modal_filters) > 0) {
          full_data <- data_reactive()
          req(full_data)
          filtered_data <- apply_filters_to_df(full_data, modal_filters)
          n_variants <- length(unique(filtered_data$variant_id))
    
          shinyWidgets::show_toast(
            title = paste(tools::toTitleCase(table_name), "Filter Applied"),
            text = paste(n_variants, "variants pass this filter. The table has been updated."),
            type = "info",
            timer = 5000,
            position = "bottom-end"
          )
        } else {
          shinyWidgets::show_toast(
            title = "Filters Cleared",
            text = "All filters for this table have been removed.",
            type = "info", timer = 5000, position = "bottom-end"
          )
        }
        removeModal()
      })
    }
    
    # handle_table_filtering("genotypes", reactive(values$table_genotypes))
    handle_table_filtering("variants", reactive(values$table_variants))
    handle_table_filtering("annotations", reactive(values$table_annotations))


    observeEvent(input$refine_plot_btn, {
      # Collect variant IDs from each filtered table.
      id_sets <- list()

      if (!is.null(values$filters_genotypes) && length(values$filters_genotypes) > 0) {
        req(values$table_genotypes)
        id_sets$geno <- apply_filters_to_df(values$table_genotypes, values$filters_genotypes)$variant_id
      }

      if (!is.null(values$filters_variants) && length(values$filters_variants) > 0) {
        req(values$table_variants)
        id_sets$vars <- apply_filters_to_df(values$table_variants, values$filters_variants)$variant_id
      }

      if (!is.null(values$filters_annotations) && length(values$filters_annotations) > 0) {
        req(values$table_annotations)
        id_sets$ann <- apply_filters_to_df(values$table_annotations, values$filters_annotations)$variant_id
      }

      if (length(id_sets) > 0) {
        # Intersect IDs to get the final variant set.
        final_ids <- Reduce(intersect, id_sets)
        values$final_selected_variants <- final_ids
        # Snapshot filters for plot rendering.
        values$plot_filters_annotations <- values$filters_annotations
        values$plot_filters_genotypes <- values$filters_genotypes
        values$plot_filters_variants <- values$filters_variants
        show_toast_success(paste(length(final_ids), "variants selected after filtering."))
      } else {
        # No filters active; show all variants.
        values$final_selected_variants <- NULL
        values$plot_filters_annotations <- NULL; values$plot_filters_genotypes <- NULL; values$plot_filters_variants <- NULL
        show_toast_success("No active filters. Showing all variants.", type = "info")
      }
    })

    observeEvent(input$clear_all_filters_btn, {
      values$filters_genotypes <- NULL; values$filters_variants <- NULL; values$filters_annotations <- NULL
      values$plot_filters_annotations <- NULL; values$plot_filters_genotypes <- NULL; values$plot_filters_variants <- NULL
      values$final_selected_variants <- NULL
      filter_states$genotypes <- NULL; filter_states$variants <- NULL; filter_states$annotations <- NULL
      show_toast_success("All filters cleared.", type = "info")
    })

    observeEvent(input$design_markers_from_region_btn, {
      # Set marker candidates from filtered variants or all variants.
      marker_candidates <- if (!is.null(values$final_selected_variants) && length(values$final_selected_variants) > 0) {
        values$final_selected_variants
      } else {
        # Fallback to all variants in the region if no filters are applied.
        req(values$table_variants)
        unique(values$table_variants$variant_id)
      }

      if (length(marker_candidates) == 0) {
        shinyWidgets::show_alert(
          title = "No Variants to Design",
          text = "There are no variants in the current filtered selection to proceed with marker design.",
          type = "warning"
        )
        return()
      }

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Preparing variants for marker design..."
      )

      c_type <- rv$conn_type
      d_path <- rv$db_path
      meta_cols <- c("chrom", "pos", "ref", "alt", "variant_type", 'major_allele', "minor_allele", "major_allele_freq", "minor_allele_freq")

      p <- future::future({
        if (c_type == "sqlite") {
          temp_con <- connect_local_db(folder_path = d_path, quiet = TRUE)
          on.exit(disconnect_local_db(temp_con, quiet = TRUE))
          fetch_genotypes_by_id(con = temp_con, variant_ids = marker_candidates, meta_data = meta_cols)
        } else {
          fetch_genotypes_by_id(connect_db_mode = 'online', variant_ids = marker_candidates, meta_data = meta_cols)
        }
      }, seed = TRUE, packages = c("panGenomeBreedr"))

      promises::then(
        p,
        onFulfilled = function(geno_data) {
          shinybusy::remove_modal_spinner()
          if (is.null(geno_data) || nrow(geno_data) == 0) {
            shinyWidgets::show_alert(title = "Error", text = "Could not retrieve genotype data for the selected variants.", type = "danger")
            return()
          }
          
          values$data_for_marker_design <- geno_data

          bslib::toggle_sidebar(id = 'db_sidebar', open = 'closed')
          update_sidebar_buttons("design_kasp_sidebar_btn")
          updateTabsetPanel(inputId = 'param_header', selected = "mark_design")
          
          updateSelectizeInput(session, inputId = "modal_marker_ID", choices = marker_candidates, selected = marker_candidates, server = TRUE)
          
          show_toast_success(paste("Loaded", length(marker_candidates), "variants for marker design."))
        },
        onRejected = function(e) {
          shinybusy::remove_modal_spinner()
          shinyWidgets::show_alert(title = "Error", text = e$message, type = "error")
        }
      )
    })

    # Design KASP Markers (from sidebar)
    observeEvent(input$design_kasp_sidebar_btn, {
      req(values$data_for_marker_design)

      bslib::toggle_sidebar(id = 'db_sidebar', open = 'closed')
      updateTabsetPanel(inputId = 'param_header', selected = "mark_design")
      updateSelectizeInput(
        session,
        inputId = "modal_marker_ID",
        choices = values$data_for_marker_design$variant_id,
        server = TRUE
      )
      update_sidebar_buttons("design_kasp_sidebar_btn")
    })

    output$filtered_variant_summary_text <- renderText({
      total_in_region <- length(unique(c(values$table_genotypes$variant_id, values$table_variants$variant_id, values$table_annotations$variant_id)))
      if (!is.null(values$final_selected_variants)) {
        paste("Visualizing", length(values$final_selected_variants), "of", total_in_region, "variants after filtering.")
      } else {
        paste("Visualizing all", total_in_region, "variants in region. No plot filters applied.")
      }
    })

    observe({
      shinyjs::toggleClass("filter_genotypes_btn", "btn-warning", !is.null(values$filters_genotypes) && length(values$filters_genotypes) > 0)
      shinyjs::toggleClass("filter_variants_btn", "btn-warning", !is.null(values$filters_variants) && length(values$filters_variants) > 0)
      shinyjs::toggleClass("filter_annotations_btn", "btn-warning", !is.null(values$filters_annotations) && length(values$filters_annotations) > 0)
    })

    # ==========================================================================
    # TABLE RENDERING LOGIC
    # ==========================================================================
    
    render_table_ui <- function(ui_id, data_slot, table_output_id) {
      output[[ui_id]] <- renderUI({
        req(values[[data_slot]])
        reactable::reactableOutput(ns(table_output_id))
      })
    }
    
    render_table_output <- function(output_id, display_reactive) {
      output[[output_id]] <- reactable::renderReactable({
        req(display_reactive())
        render_reactable(display_reactive())
      })
    }
    
    # Genotype display reactive: filters both rows (variants) and columns (samples).
    display_genotypes <- reactive({
      req(values$table_genotypes)
      data <- values$table_genotypes

      # Filter by sample metadata (columns).
      sample_ids <- values$filtered_sample_ids
      if (!is.null(sample_ids)) {
        # Preserve non-sample metadata columns.
        standard_meta_cols <- c("variant_id", "chrom", "pos", "ref", "alt", "variant_type", "major_allele", "minor_allele", "major_allele_freq", "minor_allele_freq")
        meta_cols_in_data <- intersect(standard_meta_cols, names(data))
        cols_to_keep <- c(meta_cols_in_data, intersect(sample_ids, names(data)))
        data <- data[, cols_to_keep, drop = FALSE]
      }

      # Filter by variant properties (rows).
      filters <- values$filters_genotypes
      apply_filters_to_df(data, filters)
    })

    # Generic display reactive for tables with only row filtering.
    create_display_reactive <- function(table_slot, filter_slot) {
      reactive({
        req(values[[table_slot]])
        data <- values[[table_slot]]
        filters <- values[[filter_slot]]
        apply_filters_to_df(data, filters)
      })
    }
    display_variants <- create_display_reactive("table_variants", "filters_variants")
    display_annotations <- create_display_reactive("table_annotations", "filters_annotations")
    
    render_table_ui("regional_genotypes_ui", "table_genotypes", "regional_genotypes_table"); render_table_output("regional_genotypes_table", display_genotypes)
    render_table_ui("regional_variants_ui", "table_variants", "regional_variants_table"); render_table_output("regional_variants_table", display_variants)
    render_table_ui("regional_annotations_ui", "table_annotations", "regional_annotations_table"); render_table_output("regional_annotations_table", display_annotations)

    # ==========================================================================
    # GENOTYPE TABLE FILTERING LOGIC (SPECIALIZED)
    # ==========================================================================
    observe({
      req(rv$sample_metadata)
      meta_cols <- colnames(rv$sample_metadata)
      # Exclude non-categorical columns from filter choices.
      filterable_cols <- setdiff(meta_cols, c("lib", "array_index", "lat", "lon", "latitude", "longitude"))
      updateSelectizeInput(
        session,
        "meta_filter_col",
        choices = filterable_cols,
        selected = ""
      )
    })

    observeEvent(input$filter_genotypes_btn, {
      req(values$table_genotypes)
      df <- values$table_genotypes
      
      showModal(modalDialog(
        title = "Filter Genotype Table",
        size = "xl", easyClose = TRUE,
        bslib::navset_card_tab(
          id = ns("genotype_filter_tabs"),
          bslib::nav_panel(
            "Filter Variants (Rows)",
            div(
              class = "p-3 mb-3 bg-light rounded border",
              h5("Add New Variant Filter"),
              bslib::layout_columns(
                col_widths = c(4, 5, 3),
                selectizeInput(ns("filter_col_genotypes"), "Column", choices = NULL),
                uiOutput(ns("filter_widget_genotypes")),
                div(class = "d-flex align-items-end", actionButton(ns("add_filter_genotypes"), "Add Filter", icon = icon("plus"), class = "btn-primary w-100"))
              )
            ),
            h5("Active Variant Filters"),
            uiOutput(ns("active_filters_ui_genotypes"))
          ),
          bslib::nav_panel(
            "Filter Samples (Columns)",
            selectizeInput(ns("modal_meta_filter_col"), "Filter by:", choices = NULL, multiple = TRUE, width = "100%"),
            uiOutput(ns("modal_meta_filter_values_ui"))
          )
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("clear_genotype_filters_btn"), "Clear All Filters", class = "btn-outline-danger"),
          actionButton(ns("apply_genotype_filters_btn"), "Apply & Close", icon = icon("check"), class = "btn-success")
        )
      ))

      # Update row filter choices.
      filterable_cols <- names(df)[sapply(df, function(x) !is.list(x) && !is.data.frame(x))]
      updateSelectizeInput(session, "filter_col_genotypes", choices = filterable_cols, server = TRUE)
      
      # Update column (sample metadata) filter choices.
      req(rv$sample_metadata)
      meta_cols <- colnames(rv$sample_metadata)
      filterable_meta_cols <- setdiff(meta_cols, c("lib", "array_index", "lat", "lon", "latitude", "longitude"))
      updateSelectizeInput(session, "modal_meta_filter_col", choices = filterable_meta_cols, server = TRUE)
    })

    # Row (variant) filters
    output$filter_widget_genotypes <- renderUI({
      req(input$filter_col_genotypes)
      col <- input$filter_col_genotypes; col_data <- values$table_genotypes[[col]]
      if (is.numeric(col_data)) {
        rng <- range(col_data, na.rm = TRUE)
        sliderInput(ns("filter_val_genotypes"), "Range", min = rng[1], max = rng[2], value = rng)
      } else {
        choices <- sort(unique(na.omit(as.character(col_data))))
        selectizeInput(ns("filter_val_genotypes"), "Values", choices = choices, multiple = TRUE, options = list(placeholder = "Select value(s)..."))
      }
    })
    observeEvent(input$add_filter_genotypes, {
      req(input$filter_col_genotypes, input$filter_val_genotypes)
      new_filter <- list(column = input$filter_col_genotypes, value = input$filter_val_genotypes)
      filter_states$genotypes <- c(filter_states$genotypes, list(new_filter))
    })
    output$active_filters_ui_genotypes <- renderUI({
      filters <- filter_states$genotypes
      if (is.null(filters) || length(filters) == 0) return(tags$p("No variant filters added yet.", class = "text-muted"))
      # Reuses the generic UI logic but with specific IDs.
      tagList(lapply(seq_along(filters), function(i) {
        f <- filters[[i]]; val_str <- paste(f$value, collapse = ", "); if (is.numeric(f$value)) val_str <- paste(f$value[1], "-", f$value[2])
        div(class = "d-flex justify-content-between align-items-center p-2 border-bottom",
            tags$span(tags$b(f$column), ": ", val_str),
            actionButton(ns(paste0("remove_filter_genotypes_", i)), "", icon = icon("xmark"), class = "btn-sm btn-outline-danger"))
      }))
    })
    observe({
      filters <- filter_states$genotypes
      if (!is.null(filters) && length(filters) > 0) {
        lapply(seq_along(filters), function(i) {
          observeEvent(input[[paste0("remove_filter_genotypes_", i)]], {
            current_filters <- filter_states$genotypes; current_filters[[i]] <- NULL
            filter_states$genotypes <- Filter(Negate(is.null), current_filters)
          }, ignoreInit = TRUE, once = TRUE)
        })
      }
    })

    # Column (sample) filters
    output$modal_meta_filter_values_ui <- renderUI({
      req(input$modal_meta_filter_col)
      lapply(input$modal_meta_filter_col, function(col) {
        unique_values <- sort(unique(na.omit(rv$sample_metadata[[col]])))
        selectizeInput(ns(paste0("modal_meta_vals_", col)), label = paste("Select values for:", col), choices = unique_values, multiple = TRUE, width = "100%")
      })
    })

    # Modal actions
    observeEvent(input$clear_genotype_filters_btn, {
      filter_states$genotypes <- NULL
      filter_states$sample_filters <- NULL
      updateSelectizeInput(session, "modal_meta_filter_col", selected = "")
    })

    observeEvent(input$apply_genotype_filters_btn, {
      # Apply variant filters.
      values$filters_genotypes <- filter_states$genotypes
      
      # Apply sample filters.
      if (!is.null(input$modal_meta_filter_col) && length(input$modal_meta_filter_col) > 0) {
        filters_list <- lapply(input$modal_meta_filter_col, function(col) input[[paste0("modal_meta_vals_", col)]])
        names(filters_list) <- input$modal_meta_filter_col
        filters_list <- filters_list[sapply(filters_list, function(x) !is.null(x) && length(x) > 0)]

        if (length(filters_list) > 0) {
          filtered_meta <- rv$sample_metadata
          for (col_name in names(filters_list)) {
            filtered_meta <- filtered_meta[filtered_meta[[col_name]] %in% filters_list[[col_name]], , drop = FALSE]
          }
          values$filtered_sample_ids <- filtered_meta$lib
        } else {
          values$filtered_sample_ids <- NULL
        }
      } else {
        values$filtered_sample_ids <- NULL
      }
      
      show_toast_success("Genotype table filters applied.", type = "info")
      removeModal()
    })

    observeEvent(input$clear_meta_filter, {
      values$filtered_sample_ids <- NULL
      # Clear UI.
      updateSelectizeInput(session, "meta_filter_col", selected = "")
      show_toast_success("Sample filters cleared.", type = "info")
    })

    output$filtered_sample_count <- renderText({
      if (!is.null(values$filtered_sample_ids)) {
        total_samples <- length(unique(rv$sample_metadata$lib))
        paste("Showing", length(values$filtered_sample_ids), "of", total_samples, "samples.")
      } else if (!is.null(values$table_genotypes)) {
        # Count sample columns by excluding known metadata columns.
        standard_meta_cols <- c("variant_id", "chrom", "pos", "ref", "alt", "variant_type", "major_allele", "minor_allele", "major_allele_freq", "minor_allele_freq")
        sample_cols <- setdiff(names(values$table_genotypes), standard_meta_cols)
        total_samples <- length(sample_cols)

        paste0("Showing all ", total_samples, " samples.")
      } else {
        "No samples to display."
      }
    })

    # UI and Server for Annotation Summary
    output$annotation_summary_ui <- renderUI({
      req(values$table_annotation_summary)
      tagList(
        bslib::card(
          class = "shadow-sm mb-3",
          bslib::card_header("Annotation Type Summary"),
          bslib::card_body(
            reactable::reactableOutput(ns("ann_summary_tbl"))
          )
        ),
        bslib::card(
          class = "shadow-sm mb-3",
          bslib::card_header("Impact Summary"),
          bslib::card_body(
            reactable::reactableOutput(ns("impact_summary_tbl"))
          )
        ),
        bslib::card(
          class = "shadow-sm",
          bslib::card_header("Variant Type Totals"),
          bslib::card_body(
            reactable::reactableOutput(ns("variant_totals_tbl"))
          )
        )
      )
    })

    output$ann_summary_tbl <- reactable::renderReactable({
      req(values$table_annotation_summary)
      render_reactable(values$table_annotation_summary$annotation_summary)
    })
    output$impact_summary_tbl <- reactable::renderReactable({
      req(values$table_annotation_summary)
      render_reactable(values$table_annotation_summary$impact_summary)
    })
    output$variant_totals_tbl <- reactable::renderReactable({
      req(values$table_annotation_summary)
      render_reactable(values$table_annotation_summary$variant_type_totals)
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

    # Update LD target variants when filter results change.
    observe({
      # Trigger only when LD panel is visible.
      req(
        "Linkage Disequilibrium Analysis" %in% input$explore_data_accordion,
        values$table_genotypes
      )

      # React to changes in the final, plotted variant set.
      final_variants <- values$final_selected_variants
      
      variant_choices <- values$table_genotypes$variant_id
      
      # Pre-select the variants shown in the plot.
      pre_selected_variants <- if (!is.null(final_variants)) {
        # Ensure pre-selections are valid choices.
        intersect(final_variants, variant_choices)
      } else {
        NULL
      }
      
      updateSelectizeInput(
        session,
        "target_ids",
        choices = variant_choices,
        selected = pre_selected_variants,
        server = TRUE
      )
    })

    output$ld_analysis_ui <- renderUI({
      req(values$table_genotypes)

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
      req(display_genotypes())
      geno_data <- display_genotypes()

      if (nrow(geno_data) < 2) {
        shinyWidgets::show_alert(
          title = "Not enough variants",
          text = "LD analysis requires at least two variants in the current selection.",
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

          ld_res <- calculate_LD(
            df = geno_data,
            target_variant_ids = target_variants
          )

          plot_pkg <- plot_ld_geodesic(
            ld_df = ld_res,
            query_db_geno = geno_data, # Use the filtered geno data for positions
            query_db_annot = values$table_annotations, # Use full annotations for impact lookup
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

    observeEvent(input$go_back, {
      updateTabsetPanel(session, inputId = 'param_header', selected = 'gene_cord')
      bslib::toggle_sidebar(id = 'db_sidebar', open = 'open')
      update_sidebar_buttons("show_gene_cord_btn")
    })

    # ==========================================================================
    # KASP MARKER DESIGN LOGIC
    # ==========================================================================

    kasp_des.result <- reactiveVal(NULL)
    kasp_des.plot <- reactiveVal(NULL)

    observeEvent(input$modal_run_but, {
      # Isolate target alleles for locus assay development
      active_data <- values$data_for_marker_design
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