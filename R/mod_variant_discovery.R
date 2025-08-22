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
#' @importFrom shiny NS tagList conditionalPanel textInput actionButton
#' @importFrom shiny uiOutput selectInput numericInput verbatimTextOutput
#' @importFrom shiny icon div h4 span helpText tags
#' @importFrom bslib navset_card_underline sidebar card card_header
#' @importFrom bslib card_body card_footer nav_panel navset_card_tab
#' @importFrom bslib layout_column_wrap input_switch
#'
mod_variant_discovery_ui <- function(id) {
  ns <- NS(id)

  # Reusable card styling.
  card_style <- "max-width: 350px;"
  primary_card_header <- function(title) {
    bslib::card_header(class = "bg-primary text-white", title)
  }

  # Connection Panel Components
  connection_panel <- function(ns) {
    tagList(
      div(
        id = ns("connection_panel"),
        icon("database"),
        uiOutput(ns("status_badge")),
        style = "display: flex; align-items: center; margin-bottom: 15px;"
      ),
      conditionalPanel(
        condition = paste0("output['", ns("is_connected"), "'] == false"),
        shinyFiles::shinyFilesButton(
          id = ns("connect_btn"),
          label = "Connect to Database",
          icon = icon("plug"), title = "",
          multiple = FALSE,
          style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
          `onmouseover` = "this.style.backgroundColor='#145214'",
          `onmouseout` = "this.style.backgroundColor='forestgreen'",
          buttonType = 'btn-primary'

        )
      ),
      conditionalPanel(
        condition = paste0("output['", ns("is_connected"), "'] == true"),
        actionButton(ns("disconnect_btn"), "Disconnect",
          icon = icon("unlink"),
          width = "100%",
          class = "btn-danger"
        ),
        bslib::card(
          bslib::card_header("Database Overview"),
          bslib::card_body(
            verbatimTextOutput(ns("db_info"))
          )
        )
      )
    )
  }

  # Info Tab Components
  info_tab <- function(ns) {
    conditionalPanel(
      condition = paste0("output['", ns("is_connected"), "'] == true"),
      tagList(
        bslib::accordion(
          style = "margin-bottom: 10px;",
          open = TRUE,
          bslib::accordion_panel(
            "Variant Impact Summary",
            reactable::reactableOutput(ns("table_impact_id"))
          )
        ),
        bslib::accordion(
          style = "margin-bottom: 10px;",
          open = TRUE,
          bslib::accordion_panel(
            "Variant Statistics",
            reactable::reactableOutput(ns("table_var_stats_id"))
          )
        ),
        bslib::accordion(
          style = "margin-bottom: 10px;",
          open = TRUE,
          splitLayout(
            bslib::accordion_panel(
              "Summarised SQLite Tables",
              reactable::reactableOutput(ns("sum_sqlite_id"))
            ),
            bslib::accordion_panel(
              "Variant Type Count",
              reactable::reactableOutput(ns("count_variant_typ_id"))
            )
          )
        ),
        bslib::accordion(
          style = "margin-bottom: 5px;",
          open = TRUE,
          bslib::accordion_panel(
            "Inspect SQLite Table Schema",
            selectInput(
              inputId = ns("table_name_lst"),
              label = "Table to Query",
              choices = c("variants", "annotations", "genotypes"),
              selected = "genotypes"
            ),
            reactable::reactableOutput(ns("results_lst"))
          )
        )
      )
    )
  }


  # Query Actions Tab Components
  query_actions_tab <- function(ns) {
    tagList(
      conditionalPanel(
        condition = paste0("output['", ns("is_connected"), "'] == false"),
        bslib::card(
          bslib::card_header("Variant Discovery"),
          bslib::card_body(
            "Please connect to a database to use this functionality."
          )
        )
      ),
      conditionalPanel(
        condition = paste0("output['", ns("is_connected"), "'] == true"),
        bslib::layout_column_wrap(
          width = 1 / 2,
          query_database_card(ns),
          query_action_card(ns)
        ),
        bslib::navset_card_tab(
          id = ns("nav_id"), selected = "Main Database Results ",
          bslib::nav_panel(
            title = "Main Database Results ",
            uiOutput(ns("query_db_display"))
          ),
          bslib::nav_panel(
            title = "PCVs for KASP Marker Design",
            uiOutput(ns("pcvs_kasp_marker_design_result"))
          )
        )
      )
    )
  }

  # Query Database Card Component
  query_database_card <- function(ns) {
    # Card for genotype cordinate management and database query.
    bslib::card(
      bslib::card_header(tags$strong("Set Genomic Region & Query Database")),
      tagList(
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            inputId = ns("get_cord"),
            label = "Load Genomic Region from GFF",
            icon = icon("cloud-arrow-down"),
            width = "70%"
          )
        ),
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            inputId = ns("set_cord"),
            icon = icon("pencil"),
            label = "Enter Genomic Region Manually",
            width = "70%"
          )
        )
      ),
      radioButtons(
        width = "100%",
        inputId = ns("query_database"), label = "",
        choices = c(
          "Query Database by Coordinates" = "q_entire",
          "View Annotation Summary for Region" = "q_annt"
        ),
        selected = character(0)
      ),
      uiOutput(outputId = ns("display_qd_choice")), # UI for input parameter for users based on choice
      bslib::card_footer(
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            inputId = ns("query_dbase_btn"),
            label = tags$b("Query"),
            width = "70%",
            icon = icon("database"),
            # class = "btn-info"
            style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
            `onmouseover` = "this.style.backgroundColor='#145214'",
            `onmouseout` = "this.style.backgroundColor='forestgreen'"
          )
        )
      )
    )
  }


  # Query Action Card Component
  query_action_card <- function(ns) {
    bslib::card(
      bslib::card_header(tags$strong("Filter Putative Causal Variants")),
      # Input widget for Impact levels
      selectInput(
        inputId = ns("impact_level"), label = "Impact Level",
        choices = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
        width = "100%"
      ),
      # Input widget for allele frequency
      sliderInput(
        inputId = ns("af_range"),
        label = "Allele Frequency Range",
        min = 0,
        max = 1,
        value = c(0.05, 0.9),
        step = 0.01,
        width = "100%"
      ),
      bslib::card_footer(
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            inputId = ns("get_pcv_btn"),
            label = tags$b("Extract PCVs"),
            class = "btn-warning",
            width = "70%",
            icon = icon("filter")
          )
        )
      )
    )
  }

  # Main UI Structure
  tagList(
    shinyjs::useShinyjs(),
    bslib::navset_card_underline(
      id = ns('param_header'),
      sidebar = bslib::sidebar(
        width = 350,
        title = "Database Connection",
        status = "primary",
        connection_panel(ns)
      ),
      # Display query frontend
      bslib::nav_panel(
        value = 'query_tab',
        title = tags$b("Query Actions"), icon = icon("bolt"),
        query_actions_tab(ns)
      ),
      # Database info ui
      bslib::nav_panel(
        title = tags$b("Info"), icon = icon("info-circle"),
        shinyWidgets::actionBttn(
          inputId = ns("get_info"),
          label = "Get Database Info",
          style = "float",
          size = "md",
          color = "primary"
        ),
        info_tab(ns)
      ),
      bslib::nav_panel_hidden(
        value = "mark_design",
          fluidRow(
            column(
              width = 4,
              bslib::card(
                # bslib::card_header(
                #   h4(tags$b("Input Parameters"))
                #
                #   # class = "bg-primary text-white"
                # ),
                bslib::card_body(
                  fileInput(
                    ns("modal_genome_file"),
                    label = "Genome Reference File",
                    accept = c(".fa", ".fasta", ".gz"),
                    width = "100%"
                  ),
                  selectizeInput(
                    ns("modal_marker_ID"),
                    label = "Marker ID",
                    choices = NULL,
                    options = list(placeholder = "Select variants..."),
                    multiple = TRUE,
                    width = "100%"
                  ),
                  textInput(
                    ns("modal_reg_name"),
                    label = "Region Name",
                    width = "100%",
                    placeholder = "lgs1"
                  ),
                  numericInput(
                    ns("modal_maf"),
                    width = "100%",
                    label = "Minor Allele Frequency (MAF)",
                    value = 0.05, min = 0, max = 1, step = 0.01
                  ),
                  bslib::input_switch(
                    ns("modal_draw_plot"),
                    width = "100%",
                    label = "Generate Alignment Plot",
                    value = TRUE
                  )
                ),
                bslib::card_footer(
                  actionButton(
                    width = "100%",
                    ns("modal_run_but"),
                    label = "Design Marker",
                    icon = icon("drafting-compass"),
                    # class = "btn-info"
                    style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
                    `onmouseover` = "this.style.backgroundColor='#145214'",
                    `onmouseout` = "this.style.backgroundColor='forestgreen'"
                  )
                )
              )
            ),
            column(
              width = 8,
              bslib::accordion(
                style = "margin-bottom: 70px;",
                id = ns("results_accordion"),
                width = "100%",
                open = TRUE,
                bslib::accordion_panel(
                  "KASP Marker Data & Sequence Alignment Table",
                  DT::DTOutput(ns("kasp_table")),
                    fluidRow(
                      column(width = 3, selectInput(
                        inputId = ns("exten"),
                        label = "Download file as?",
                        choices = c(".csv", ".xlsx"),
                        selected = ".xlsx",
                        multiple = FALSE,
                      )),
                      column(
                        width = 4,
                        textInput(
                          inputId = ns("file_name"),
                          label = "Enter File Prefix",
                          value = "Kasp M_D for Intertek"
                        )
                      )
                    ),
                    downloadButton(ns("download_table"),
                      label = "Download File",
                      class = "btn-success",
                      icon = icon("download")
                    )
                )
                ),
                uiOutput(ns("plot_container"))

            )
          ),
            actionButton(
              inputId = ns("go_back"),
              width = '10%',
              label   = "Back",
              icon    = icon("arrow-left"),
              class   = "btn-secondary"
            )
      )
    )
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

    # Reactive values to store application state
    rv <- reactiveValues(
      conn = NULL, # Connection object
      connected = FALSE,
      tables = NULL,
      status_message = "Not connected",
      selected_table = NULL,
      variant_impact = NULL,
      db_path = NULL,
      sqlite_summary = NULL,
      variant_count = NULL,
      variant_stats = NULL,
      lst_tbl_column = NULL
    )

    # Status badge UI
    output$status_badge <- renderUI({
      status_class <- if (rv$connected) "status-connected" else "status-disconnected"
      status_icon <- if (rv$connected) "check-circle" else "times-circle"
      status_text <- if (rv$connected) " Connected" else " Disconnected"

      span(
        class = status_class, style = "margin-left: 10px;",
        icon(status_icon), status_text
      )
    })


    # DATABASE CONNECTION MANAGEMENT
    #--------------------------------------------
    volumes <- shinyFiles::getVolumes() # extract all paths
    # Setup file chooser
    shinyFiles::shinyFileChoose(input, "connect_btn", roots = volumes, session = session)

    # Connect to database
    observeEvent(input$connect_btn, {
      file_selected <- shinyFiles::parseFilePaths(volumes, input$connect_btn)

      if (nrow(file_selected) > 0) {

        shinybusy::show_modal_spinner(
          spin = "fading-circle",
          color = "#0dc5c1",
          text = "Connecting to Database... Please wait."
        )

        tryCatch({
          # Extract file path
          db_path <- as.character(file_selected$datapath)

        if (!file.exists(db_path)) {
          shinyWidgets::show_toast("Error: Database file not found", type = "error")
          return()
        }

        if (!is.null(rv$conn)) {
          dbDisconnect(rv$conn)
          rv$conn <- NULL
        }

        rv$conn <- dbConnect(RSQLite::SQLite(), db_path)
        rv$connected <- TRUE
        rv$tables <- list_sqlite_tables(db_path)
        rv$status_message <- paste("Connected with", length(rv$tables), "tables")
        rv$db_path <- db_path

        shinyWidgets::show_alert(
          title = "Success!",
          text = "Database connected successfully",
          type = "success",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Failed!",
          text = "Unable to connect to database",
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
      }
    })


    # Disconnect from database
    observeEvent(input$disconnect_btn, {
      if (!is.null(rv$conn)) {
        dbDisconnect(rv$conn)
        rv$conn <- NULL

        # Reset state
        rv$connected <- FALSE
        rv$tables <- NULL
        rv$status_message <- "Disconnected"
        rv$variant_impact <- NULL
        rv$selected_table <- NULL
        rv$sqlite_summary <- NULL
        rv$variant_count <- NULL
        rv$variant_stats <- NULL
        rv$lst_tbl_column <- NULL
        rv$db_path <- NULL

        # Reset query results
        values$result <- NULL
        values$query_db_val <- NULL
        values$query_ann_react <- NULL
        values$query_geno_react <- NULL

        shinyWidgets::show_toast(
          title = "Disconnected from Database",
          type = "warning",
          timer = 5000
        )
      }
    })

    # Clean up connection when app closes
    session$onSessionEnded(function() {
      isolate({
        if (!is.null(rv$conn)) {
          dbDisconnect(rv$conn)
          rv$conn <- NULL
        }
      })
    })

    # Connection status for conditional panels
    output$is_connected <- reactive(rv$connected)
    outputOptions(output, "is_connected", suspendWhenHidden = FALSE)

    # Database information output
    output$db_info <- renderPrint({
      req(rv$connected)
      cat("Path: ", rv$db_path, "\n")
      cat("Tables: ", length(rv$tables), "\n")
      cat("Available tables:\n")
      cat(paste("- ", rv$tables, collapse = "\n"))
    })


    # Info section once database is loaded #
    #--------------------------------------

    # Users decide to get info about table then it renders.
    observeEvent(input$get_info, {
      req(rv$db_path)
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Retrieving database information"
      )

      tryCatch({
        # Precompute database information
        rv$variant_impact <- variant_impact_summary(db_path = rv$db_path)
        rv$sqlite_summary <- summarize_sqlite_tables(db_path = rv$db_path)
        rv$variant_count <- count_variant_types(db_path = rv$db_path)
        rv$variant_stats <- variant_stats(db_path = rv$db_path)
        rv$tables <- list_sqlite_tables(rv$db_path)
        rv$lst_tbl_column <- list_table_columns(
          db_path = rv$db_path,
          table_name = input$table_name_lst
        )
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Failed!",
          text = e$message,
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })

    # Table renderers
    # Helper function for rendering tables
    render_reactable <- function(data, theme = NULL) {
      # Default theme if none provided
      if (is.null(theme)) {
        theme <- reactable::reactableTheme(
          backgroundColor = "hsl(0, 0%, 100%)",
          borderColor = "hsl(0, 0%, 89%)",
          stripedColor = "hsl(0, 0%, 97%)",
          highlightColor = "hsl(0, 0%, 96%)",
          cellPadding = "8px 12px",
          style = list(
            fontFamily = "-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif"
          ),
          searchInputStyle = list(
            width = "100%",
            padding = "8px 12px"
          ),
          headerStyle = list(
            borderWidth = "1px",
            padding = "8px 12px"
          )
        )
      }

      reactable::reactable(
        data,
        defaultColDef = reactable::colDef(
          headerClass = "header",
          align = "left",
          minWidth = 100,
          width = 200
        ),
        pagination = TRUE,
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(10, 25, 50, 100),
        striped = TRUE,
        highlight = TRUE,
        bordered = TRUE,
        compact = TRUE,
        wrap = TRUE,
        theme = theme,
        style = list(
          maxWidth = "100%"
        )
      )
    }

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




    # Query section, UI rendering.
    #------------------------------
    # Reactive values for query results
    values <- reactiveValues(
      result = NULL, # Store results for genotype coordinates
      query_db_val = NULL, # Store results for queried tables within Database
      query_ann_react = NULL, # store annotation summary
      query_geno_react = NULL, # store genotype matrix of PCVs
      last_action = NULL # help make feedback robust for dbquery actions
    )

    # Helper function for success toasts
    show_toast_success <- function(text,
                                   type = "success",
                                   timer = 3000) {
      shinyWidgets::show_toast(
        title = "Success",
        text = text,
        timer = timer,
        position = "bottom-end",
        type = type
      )
    }

    # Helper function for annotation summary ui
    annotation_summary_results_ui <- function(ns) {
      bslib::card(
        bslib::card_header("Annotation Statistics"),
        bslib::card(
          bslib::card_header("Annotation Summary", class = "bg-primary text-light"),
          reactable::reactableOutput(ns("ann_summary_tbl"))
        ),
        bslib::card(
          bslib::card_header("Impact Summary", class = "bg-primary text-light"),
          reactable::reactableOutput(ns("impact_summary_tbl"))
        ),
        bslib::card(
          bslib::card_header("Variant Type Totals", class = "bg-primary text-light"),
          reactable::reactableOutput(ns("variant_totals_tbl"))
        ),
        class = "mt-3"
      )
    }



    # Gene coordinates modal
    observeEvent(input$get_cord, {
      showModal(
        modalDialog(
          title = div("Get Genotype Co-ordinates", style = "text-align: center;"),
          easyClose = FALSE,
          footer = tagList(
            actionButton(
              inputId = ns("dismiss_modal1"),
              label = "Cancel",
              class = "btn-danger",
              icon = icon("times")
            )
          ),
          bslib::card(
            # bslib::card_header("Gene Parameters"),
            textInput(
              inputId = ns("gene_name"), label = "Gene Name (Sobic ID)",
              value = "Sobic.005G213600",
              placeholder = "Enter Sobic ID (e.g., Sobic.005G213600)", width = "100%"
            ),
            radioButtons(
              inputId = ns("input_method"), label = "GFF File Source",
              choices = c("URL" = "url", "Upload File" = "file"),
              selected = "url"
            ),
            conditionalPanel(
              condition = paste0("input['", ns("input_method"), "'] === 'url'"),
              textInput(
                inputId = ns("gff_url"), label = "GFF File URL", width = "100%",
                value = "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz",
                placeholder = "Enter URL to GFF file"
              )
            ),
            conditionalPanel(
              condition = paste0("input['", ns("input_method"), "'] === 'file'"),
              fileInput(
                inputId = ns("gff_file"), label = "Upload GFF File",
                accept = c(".gff3", ".gff", ".gff3.gz", ".gff.gz"), width = "100%"
              )
            ),
            bslib::card_footer(
              actionButton(
                inputId = ns("submit"),
                label = "Get Coordinates",
                width = "100%",
                # class = "btn-info",
                style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
                `onmouseover` = "this.style.backgroundColor='#145214'",
                `onmouseout` = "this.style.backgroundColor='forestgreen'",
                icon = icon("search")
              )
            )
          )
        )
      )
    })

    # Remove modal if cancel is clicked.
    observeEvent(input$dismiss_modal1, {
      removeModal()
    })

    # Get gene coordinates
    observeEvent(input$submit, {
      removeModal()
      req(input$gene_name)

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Getting gene coordinates... Please wait."
      )

      tryCatch({
        # Assigning gff input choice correctly.
        if (input$input_method == "url") {
          req(input$gff_url)
          gff_path <- input$gff_url
        } else if (input$input_method == "file") {
          req(input$gff_file)

          gff_path <- input$gff_file$datapath
        }

        values$result <- gene_coord_gff(input$gene_name, gff_path)

        shinyWidgets::show_alert(
          title = "Found Gene Co-ordinates",
          # Informative message to confirm that cordinates have been found
          text = sprintf(
            "Chromosome: %s | Start: %d | End: %d",
            values$result$chrom, values$result$start, values$result$end
          ),
          type = "success",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Failed!",
          text = e$message,
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })



    # If user selects to manually input the cordinates.
    observeEvent(input$set_cord, {
      showModal(
        modalDialog(
          title = div(
            style = "text-align: center;",
            "Set Genotype Co-ordinates Manually"
          ),
          easyClose = FALSE,
          footer = tagList(
            actionButton(
              inputId = ns("dismiss_modal"),
              label = "Cancel",
              class = "btn-danger",
              icon = icon("times")
            )
          ),
          tagList(
            bslib::card(
              textInput(
                inputId = ns("chrom"),
                label = "Chromosome",
                value = NULL,
                placeholder = 'Chr05',
                width = "100%"
              ), # chromosome name input
              numericInput(
                inputId = ns("start"),
                label = "Start Position",
                value = NULL,
                width = "100%"
              ), # chrom start input
              numericInput(
                inputId = ns("end"),
                label = "End Position",
                value = NULL,
                width = "100%"
              ), # chrom end input
              # footer
              bslib::card_footer(
                actionButton(
                  inputId = ns("set_genocod_btn"),
                  label = "Submit",
                  width = "100%",
                  # class = "btn-info",
                  style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
                  `onmouseover` = "this.style.backgroundColor='#145214'",
                  `onmouseout` = "this.style.backgroundColor='forestgreen'",
                  icon = icon("edit")
                )
              )
            )
          )
        )
      )
    })

    # Close modal is close is clicked.
    observeEvent(input$dismiss_modal, {
      removeModal()
    })

    # Close modal after manual setting
    observeEvent(input$set_genocod_btn, {
      req(input$chrom, input$start, input$end)
      removeModal()

      values$result <- list(
        chrom = input$chrom,
        start = input$start,
        end = input$end
      )

      shinyWidgets::show_alert(
        title = "Gene Coordinates Set",
        text = sprintf(
          "Chromosome: %s | Start: %d | End: %d",
          input$chrom, input$start, input$end
        ),
        type = "success",
        showCloseButton = TRUE,
        timer = 5000
      )
    })



    # Render UI dynamically based on choice of query
    observe({
      if (!is.null(input$query_database) && input$query_database == "q_entire") {
        output$display_qd_choice <- renderUI({
          tagList(
            selectInput(
              inputId = ns("table_name"),
              label = "Table to Query",
              choices = c("variants", "annotations", "genotypes"),
              selected = "annotations", width = "100%"
            ),
            uiOutput(outputId = ns("gene_name_id"))
          ) # UI renders when table == annotations
        })
      } else if (!is.null(input$query_database) && input$query_database == "q_annt") {
        output$display_qd_choice <- renderUI({
          tagList(
            selectInput(
              inputId = ns("table_name_a"),
              label = "Table name for Annotation",
              choices = c("variants", "annotations", "genotypes"),
              selected = "annotations", width = "100%"
            ),
            selectInput(
              inputId = ns("table_name_v"),
              label = "Table name for Variants",
              choices = c("variants", "annotations", "genotypes"),
              selected = "variants", width = "100%"
            )
          )
        })
      } else if (is.null(input$query_database)) {
        output$display_qd_choice <- renderUI({
          NULL
        })
      }
    })

    # Gene name input for annotations
    output$gene_name_id <- renderUI({
      if (input$table_name == "annotations") {
        textInput(
          inputId = ns("query_gene_name"),
          label = "Gene Name",
          value = NULL,
          width = "100%"
        )
      }
    })


    # Combined database query and annotation summary
    observeEvent(input$query_dbase_btn, {
      values$last_action <- NULL # Null last action

      req(rv$db_path, values$result)

      # Show loading spinner for both operations
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Querying Database... Please wait."
      )


      # First: Database query
      tryCatch({
        if (!is.null(input$table_name) && input$table_name != "") {
          values$query_db_val <- query_db(
            db_path = rv$db_path,
            table_name = input$table_name,
            chrom = values$result$chrom,
            start = values$result$start,
            end = values$result$end,
            gene_name = if (input$query_gene_name == "") NULL else input$query_gene_name
          )

          values$last_action <- "main_db"
        }

        # Second: Annotation summary
        if (!is.null(input$table_name_v) && !is.null(input$table_name_a)) {
          updateTabsetPanel(session, "nav_id", selected = "Main Database Result")

          values$query_ann_react <- query_ann_summary(
            db_path = rv$db_path,
            variants_table = input$table_name_v,
            annotations_table = input$table_name_a,
            chrom = values$result$chrom,
            start = values$result$start,
            end = values$result$end
          )

          values$last_action <- "annotation_summ"
        }
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Failed!",
          text = e$message,
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinyjs::delay(ms = 1000, {
          shinybusy::remove_modal_spinner()
          if (values$last_action == "annotation_summ") {
            show_toast_success("Variant Filtered By Annotation Summary")
          } else if (values$last_action == "main_db") {
            show_toast_success(text = paste("Queried by", input$table_name))
          }
        })
      })
    })


    # Insert results to table
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


    # Query results display # superrendered UI
    observe({
      req(input$query_database) # choice radio button
      if (!is.null(input$query_database) && input$query_database == "q_entire") {
        output$query_db_display <- renderUI({
          reactable::renderReactable({
            req(values$query_db_val)
            render_reactable(values$query_db_val)
          })
        })
      } else if (!is.null(input$query_database) && input$query_database == "q_annt") {
        output$query_db_display <- renderUI({
          annotation_summary_results_ui(ns) # display cards showing  list of annotation summary
        })
      }
    })





    # Get Putative Causal variants section
    #----------------
    genotype_results_ui <- function(ns) {
      bslib::card(
        # bslib::card_header("Genotype Data"),
        div(style = "overflow-x: auto;", reactable::reactableOutput(ns("genotype_results_tbl"))),
        class = "mt-3",
        bslib::card_footer(
          textInput(
            inputId = ns("File_name"),
            label = "Enter File Name",
            value = "Chrom 05",
            width = "30%"
          ),
          div(
            style = "display: flex; justify-content: space-between;
            align-items: center; gap: 10px; width: 100%; flex-grow: 1;",
            shinyWidgets::downloadBttn(
              color = "primary",
              style = "unite",
              icon = icon("upload"),
              outputId = ns("download_excel"),
              label = "Export Genotype Matrix as .xlsx"
            ),
            shinyWidgets::actionBttn(
              style = "unite",
              color = "success",
              inputId = ns("push_1"),
              label = "Design KASP Markers for PCVs",
              icon = icon("dna")
            )
          )
        )
      )
    }


    # Query handlers
    # Reactive for queried by impact
    query_by_impact_result <- reactive({
      req(rv$db_path, input$impact_level, values$result)
      query_by_impact(
        db_path = rv$db_path,
        impact_level = input$impact_level,
        chrom = values$result$chrom,
        start = values$result$start,
        end = values$result$end
      )
    })

    # Hold genotypes for impact queries
    hold_genotypes_impact <- reactive({
      req(query_by_impact_result(), rv$db_path)
      query_genotypes(
        db_path = rv$db_path,
        variant_ids = query_by_impact_result()$variant_id,
        meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
      )
    })

    query_by_alf_result <- reactive({
      req(input$af_range, hold_genotypes_impact())

      filter_by_af(
        gt = hold_genotypes_impact(),
        min_af = input$af_range[1],
        max_af = input$af_range[2]
      )
    })


    # Filter by out successful PCVs for marker design
    observeEvent(input$get_pcv_btn, {
      values$query_geno_react <- NULL

      req(query_by_alf_result(), rv$db_path)

      updateTabsetPanel(session, inputId = "nav_id", selected = "PCVs for KASP Marker Design")

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Getting Putative Causal Variants... Please wait."
      )

      tryCatch(
        {
          values$query_geno_react <- query_genotypes(
            db_path = rv$db_path,
            variant_ids = query_by_alf_result()$variant_id,
            meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
          )
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Error",
            text = "No Putative Causal Variants Found!",
            type = "error",
            showCloseButton = TRUE,
            timer = 8000
          )
        }, finally = {
          shinyjs::delay(1000, {
            shinybusy::remove_modal_spinner()
            if (!is.null(values$query_geno_react)) {
              show_toast_success(text = paste("Found", nrow(values$query_geno_react), "Putative Causal Variants"))
            }
          })
        }
      )

      # Display result for successful pcvs
      output$pcvs_kasp_marker_design_result <- renderUI({
        genotype_results_ui(ns)
      })
    })


    # Result genotype matrix
    output$genotype_results_tbl <- reactable::renderReactable({
      req(values$query_geno_react)
      render_reactable(values$query_geno_react)
    })



    # Download handler
    output$download_excel <- downloadHandler(
      filename = function() {
        paste("dbquery", input$File_name, "genotypes", ".xlsx", sep = "_")
      },
      content = function(file) {
        writexl::write_xlsx(values$query_geno_react, path = file)
      }
    )

    # Design Kasp marker within variant discovery.
    #----------------------------------------------
    # Show modal window when user decides to design kasp Marker.
    observeEvent(input$push_1, {
      req(values$query_geno_react)

     updateTabsetPanel(inputId = 'param_header',selected = "mark_design")
      updateSelectizeInput(
        session,
        inputId = "modal_marker_ID",
        choices = values$query_geno_react$variant_id,
        server = TRUE
      )

    })

    # When back button is clicked, go back to query panel
    observeEvent(input$go_back,{
      updateTabsetPanel(inputId = 'param_header',selected = 'query_tab')
    })

    # Server Side of kasp marker design
    kasp_des.result <- reactiveVal(NULL) # store kasp dataframe
    kasp_des.plot <- reactiveVal(NULL) # store kasp plots

    observeEvent(input$modal_run_but, {
      req(
        values$query_geno_react, input$modal_genome_file$datapath,
        input$modal_maf, input$modal_marker_ID
      )

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Designing KASP Marker... Please wait."
      )


      list_markers <- list() # list object for markers
      list_plots <- list() # list object for plots

      tryCatch(
        {
          # Get column names
          gt_cols <- colnames(values$query_geno_react)

          # Get exact column names
          id_col <- gt_cols[grep("id", gt_cols, ignore.case = TRUE)[1]]
          chrom_col <- gt_cols[grep("chro", gt_cols, ignore.case = TRUE)[1]]
          pos_col <- gt_cols[grep("pos", gt_cols, ignore.case = TRUE)[1]]
          ref_col <- gt_cols[grep("ref", gt_cols, ignore.case = TRUE)[1]]
          alt_col <- gt_cols[grep("alt", gt_cols, ignore.case = TRUE)[1]]

          # Hard coded genotype start.
          geno_start <- 7

          # Get unique chromosomes
          unique_chr <- unique(values$query_geno_react[[chrom_col]])

          for (marker in input$modal_marker_ID) {
            # Run KASP design
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


          kasp_des.result(data.table::rbindlist(list_markers)) # marker dataframes
          kasp_des.plot(list_plots) # plots

          # Success alert for VCF processing
          shinyWidgets::show_alert(
            title = "Success!",
            text = paste(
              length(list_markers), "KASP marker(s) and",
              length(list_plots), "plot(s) designed successfully"
            ),
            type = "success",
            showCloseButton = TRUE,
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
        },finally = {
          shinybusy::remove_modal_spinner()
        }
      )
    })
    # Update drop down.
    observe({
      req(kasp_des.result(), kasp_des.plot())
      updateSelectizeInput(session,
        inputId = "done_markers",
        choices = names(kasp_des.result())
      )
      updateSelectizeInput(session,
        inputId = "plot_choice",
        choices = names(kasp_des.plot())
      )
    })


    output$kasp_table <- DT::renderDT({
      req(kasp_des.result())
      DT::datatable(
        data = kasp_des.result(),
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          serverSide = TRUE
        ),
        escape = TRUE
      )
    })


    # Download Kasp marker result(excel or csv.)
    output$download_table <- downloadHandler(
      filename = function() {
        # Clean proposed user file name and append underscores
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

    # Plot container UI - show if user selects true
    observeEvent(input$modal_draw_plot, {
      if (input$modal_draw_plot == TRUE) {
        output$plot_container <- renderUI({
          bslib::accordion(
            bslib::accordion_panel(
              "KASP Sequence Alignment Plot",
              selectizeInput(
                inputId = ns("plot_choice"),
                label = "Select Marker ID",
                width = "45%",
                choices = NULL
              ), # drop down for plots
              plotOutput(ns("plot")),
              downloadButton(ns("download_plot"),
                label = "Download Plot (pdf)",
                class = "btn-success", icon = icon("download")
              )
            )
          )
        })
      } else {
        output$plot_container <- renderUI({
          NULL
        })
      }
    })


    # Render Plot
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

        # Start PDF
        grDevices::pdf(file, width = 24, height = 9, onefile = TRUE)

        # Print plots
        if (is.list(plots)) lapply(plots, print) else print(plots)

        grDevices::dev.off() # Close PDF
      }
    )
  })
}

## To be copied in the UI
# mod_variant_discovery_ui("variant_discovery_1")

## To be copied in the server
# mod_variant_discovery_server("variant_discovery_1")
