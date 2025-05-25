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
#' @importFrom DT DTOutput
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
        textInput(ns("db_path"), "Database Path:",
          placeholder = "Enter the full path to your SQLite database"
        ),
        actionButton(ns("connect_btn"), "Connect to Database",
          icon = icon("plug"),
          width = "100%",
          class = "btn-primary"
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
      bslib::layout_column_wrap(
        width = "100%",
        bslib::card(
          primary_card_header("Variant Impact Summary"),
          bslib::card_body(DT::DTOutput(ns("table_impact_id")))
        ),
        bslib::card(
          primary_card_header("Variant Statistics"),
          bslib::card_body(DT::DTOutput(ns("table_var_stats_id")))
        )
      ),
      bslib::layout_column_wrap(
        width = 1 / 2,
        bslib::card(
          primary_card_header("Summarised SQLite Tables"),
          bslib::card_body(DT::DTOutput(ns("sum_sqlite_id")))
        ),
        bslib::card(
          primary_card_header("Variant Type Count"),
          bslib::card_body(DT::DTOutput(ns("count_variant_typ_id")))
        )
      ),
      bslib::layout_column_wrap(
        width = "100%",
        bslib::card(
          primary_card_header("Inspect SQlite Table Schema"),
          bslib::card_body(
            selectInput(
              inputId = ns("table_name_lst"),
              label = "Table to Query",
              choices = c("variants", "annotations", "genotypes"),
              selected = "genotypes"
            ),
            DT::DTOutput(ns("results_lst"))
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
          id = ns("nav_id"), selected = "Main Database Query Result",
          bslib::nav_panel(
            title = "Main Database Query Result",
            uiOutput(ns("query_db_display"))
          ),
          bslib::nav_panel(
            title = "Results Showing PCVs for KASP Marker Design",
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
      bslib::card_header(tags$strong("Genotype Coordinates & Database Query")),
      tagList(
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            inputId = ns("get_cord"), label = "Get Genotype Co-ordinates from GFF",
            width = "70%"
          )
        ),
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            inputId = ns("set_cord"), label = "Set Genotype Co-ordinates Manually",
            width = "70%"
          )
        )
      ),
      radioButtons(
        width = "100%",
        inputId = ns("query_database"), label = "",
        choices = c(
          "Query Database Using Coordinates" = "q_entire",
          "Annotation Summary for Defined Region" = "q_annt"
        ),
        selected = character(0)
      ),
      uiOutput(outputId = ns("display_qd_choice")), # UI for input parameter for users based on choice
      bslib::card_footer(
        div(
          style = "display: flex; justify-content: center;",
          actionButton(
            inputId = ns("query_dbase_btn"), label = "Query Database",
            width = "70%", icon = icon("database"), class = "btn-primary"
          )
        )
      )
    )
  }


  # Query Action Card Component
  query_action_card <- function(ns) {
    bslib::card(
      bslib::card_header(tags$strong("Get Putative Causal Variants")),
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
            inputId = ns("get_pcv_btn"), label = "Extract PCVs",
            class = "btn-primary",
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
      sidebar = bslib::sidebar(
        width = 350,
        title = "Database Connection",
        status = "primary",
        connection_panel(ns)
      ),
      bslib::nav_panel(
        title = tags$b("Query Actions"), icon = icon("bolt"),
        query_actions_tab(ns)
      ), # Display query frontend
      bslib::nav_panel(
        title = tags$b("Info"), icon = icon("info-circle"),
        shinyWidgets::actionBttn(
          inputId = ns("get_info"),
          label = "Get Database Info",
          style = "float",
          size = "md",
          color = "primary"
        ),
        info_tab(ns) # Display info
      )
    )
  )
}

#' Variant Discovery Server Function
#'
#' @description Server logic for variant discovery module.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom shinyalert shinyalert
#' @importFrom RSQLite dbConnect dbDisconnect
#' @importFrom writexl write_xlsx
#' @importFrom shinyWidgets show_alert show_toast
#' @importFrom DT DTOutput
#' @importFrom shiny showModal modalDialog
#'
#' @noRd
mod_variant_discovery_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Conditional logic
    required_pkgs <- c("writexl", "shinyalert", "readxl", "shinyWidgets", "shinybusy")

    for (pkg in required_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("The '%s' package is required for this feature
                   Please install it using install.packages(%s).", pkg, pkg),
          call. = FALSE
        )
      }
    }

    # Connection object
    conn <- NULL

    # Reactive values to store application state
    rv <- reactiveValues(
      connected = FALSE,
      tables = NULL,
      status_message = "Not connected",
      selected_table = NULL,
      variant_impact = NULL,
      db_path = NULL,
      sqlite_summary = NULL,
      variant_count = NULL,
      variant_stats = NULL
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

    # Connect to database
    observeEvent(input$connect_btn, {
      req(input$db_path)

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Connecting to Database... Please wait."
      )

      tryCatch({
        if (!file.exists(input$db_path)) {
          shinyWidgets::show_toast("Error: Database file not found", type = "error")
          return()
        }

        # Close existing connection if any
        if (!is.null(conn)) {
          dbDisconnect(conn)
          conn <- NULL
        }

        # Connect to database
        conn <<- dbConnect(RSQLite::SQLite(), input$db_path)

        # Update state
        rv$connected <- TRUE
        rv$tables <- list_sqlite_tables(input$db_path)
        rv$status_message <- paste("Connected with", length(rv$tables), "tables")
        rv$db_path <- input$db_path

        shinyjs::delay(100, {
          show_alert(
            title = "Success!",
            text = "Database connected successfully",
            type = "success",
            showCloseButton = TRUE,
            timer = 5000
          )
        })
      }, error = function(e) {
        show_alert(
          title = "Failed!",
          text = "Unable to connect to database",
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })

    # Disconnect from database
    observeEvent(input$disconnect_btn, {
      if (!is.null(conn)) {
        dbDisconnect(conn)
        conn <<- NULL

        # Reset state
        rv$connected <- FALSE
        rv$tables <- NULL
        rv$status_message <- "Disconnected"
        rv$variant_impact <- NULL
        rv$selected_table <- NULL

        shinyWidgets::show_toast(
          title = "Disconnected from Database",
          type = "warning",
          timer = 5000
        )
      }
    })

    # Clean up connection when app closes
    onSessionEnded(function() {
      if (!is.null(conn)) {
        dbDisconnect(conn)
      }
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
      req(input$db_path)
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Retrieving database information… This may take some time"
      )

      tryCatch({
        # Precompute database information
        rv$variant_impact <- variant_impact_summary(db_path = input$db_path)
        rv$sqlite_summary <- summarize_sqlite_tables(db_path = input$db_path)
        rv$variant_count <- count_variant_types(db_path = input$db_path)
        rv$variant_stats <- variant_stats(db_path = input$db_path)
        rv$tables <- list_sqlite_tables(input$db_path)
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })

    # Table renderers
    # Helper function for rendering tables
    render_dt <- function(data) {
      DT::datatable(data, options = list(scrollX = TRUE))
    }

    output$table_impact_id <- DT::renderDT({
      req(rv$variant_impact)
      render_dt(rv$variant_impact)
    })

    output$table_var_stats_id <- DT::renderDT({
      req(rv$variant_stats)
      render_dt(rv$variant_stats)
    })

    output$sum_sqlite_id <- DT::renderDT({
      req(rv$sqlite_summary)
      render_dt(rv$sqlite_summary)
    })

    output$count_variant_typ_id <- DT::renderDT({
      req(rv$variant_count)
      render_dt(rv$variant_count)
    })

    # List table columns in database
    lst_tbl_column <- reactive({
      req(input$table_name_lst, input$db_path)
      list_table_columns(
        db_path = input$db_path,
        table_name = input$table_name_lst
      )
    })

    output$results_lst <- DT::renderDT({
      req(lst_tbl_column())
      render_dt(lst_tbl_column())
    })




    # Query section, UI rendering.
    #------------------------------
    # Reactive values for query results
    values <- reactiveValues(
      # status = "",
      result = NULL, # Store results for genotype coordinates
      query_db_val = NULL, # Store results for queried tables within Database
      query_ann_react = NULL, # store annotation summary
      query_geno_react = NULL # store genotype matrix of PCVs
    )

    # Helper function for success toasts
    show_toast_success <- function(text) {
      shinyWidgets::show_toast(
        title = "Success",
        text = text,
        timer = 5000,
        position = "bottom-end",
        type = "success"
      )
    }
    # Helper function for annotation summary ui
    annotation_summary_results_ui <- function(ns) {
      bslib::card(
        bslib::card_header("Annotation Statistics"),
        bslib::card(
          bslib::card_header("Annotation Summary", class = "bg-primary text-light"),
          DT::DTOutput(ns("ann_summary_tbl"))
        ),
        bslib::card(
          bslib::card_header("Impact Summary", class = "bg-primary text-light"),
          DT::DTOutput(ns("impact_summary_tbl"))
        ),
        bslib::card(
          bslib::card_header("Variant Type Totals", class = "bg-primary text-light"),
          DT::DTOutput(ns("variant_totals_tbl"))
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
                inputId = ns("submit"), label = "Get Coordinates",
                class = "btn-primary w-100",
                icon = icon("search")
              )
            )
          )
        )
      )
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

        show_alert(
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
        show_alert(
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
          title = div(style = "text-align: center;", "Set Genotype Co-ordinates Manually"), easyClose = FALSE,
          tagList(
            bslib::card(
              textInput(
                inputId = ns("chrom"),
                label = "Chromosome",
                value = NULL,
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
                  class = "btn-primary w-100",
                  icon = icon("edit")
                )
              )
            )
          )
        )
      )
    })

    # Close modal after manual setting
    observeEvent(input$set_genocod_btn, {
      req(input$chrom, input$start, input$end)
      removeModal()

      show_alert(
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



    # Query database
    observeEvent(input$query_dbase_btn, {
      req(input$db_path, input$table_name)

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Querying Database... Please wait."
      )

      tryCatch({
        values$query_db_val <- query_db(
          db_path = input$db_path,
          table_name = input$table_name,
          chrom = if (is.null(input$chrom)) values$result$chrom else input$chrom,
          start = if (is.null(input$start)) values$result$start else input$start,
          end = if (is.null(input$end)) values$result$end else input$end,
          gene_name = if (input$query_gene_name == "") NULL else input$query_gene_name
        )

        show_toast_success(text = paste("Queried by", input$table_name))
      }, error = function(e) {
        show_alert(
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

    # Annotation summary.
    observeEvent(input$query_dbase_btn, {
      updateTabsetPanel(session, "nav_id", selected = "Main Database Query Result")
      req(input$db_path, input$table_name_v, input$table_name_a)

      values$query_ann_react <- query_ann_summary(
        db_path = input$db_path,
        variants_table = input$table_name_v,
        annotations_table = input$table_name_a,
        chrom = if (is.null(input$chrom)) values$result$chrom else input$chrom,
        start = if (is.null(input$start)) values$result$start else input$start,
        end = if (is.null(input$end)) values$result$end else input$end
      )

      show_toast_success("Variant Filtered By Annotation Summary")
    })

    # Insert results to table
    output$ann_summary_tbl <- DT::renderDT({
      req(values$query_ann_react)
      render_dt(values$query_ann_react$annotation_summary)
    })

    output$impact_summary_tbl <- DT::renderDT({
      req(values$query_ann_react)
      render_dt(values$query_ann_react$impact_summary)
    })

    output$variant_totals_tbl <- DT::renderDT({
      req(values$query_ann_react)
      render_dt(values$query_ann_react$variant_type_totals)
    })


    # Query results display # superrendered UI
    observe({
      req(input$query_database) # choice radio button
      if (!is.null(input$query_database) && input$query_database == "q_entire") {
        output$query_db_display <- renderUI({
          DT::renderDT({
            req(values$query_db_val)
            render_dt(values$query_db_val)
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
        div(style = "overflow-x: auto;", DT::DTOutput(ns("genotype_results_tbl"))),
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
      req(input$db_path, input$impact_level)
      query_by_impact(
        db_path = input$db_path,
        impact_level = input$impact_level,
        chrom = if (is.null(input$chrom)) values$result$chrom else input$chrom,
        start = if (is.null(input$start)) values$result$start else input$start,
        end = if (is.null(input$end)) values$result$end else input$end
      )
    })

    # Hold genotypes for impact queries
    hold_genotypes_impact <- reactive({
      req(query_by_impact_result(), input$db_path)
      query_genotypes(
        db_path = input$db_path,
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
      updateTabsetPanel(session, inputId = "nav_id", selected = "Results Showing PCVs for KASP Marker Design")

      req(query_by_alf_result(), input$db_path)

      tryCatch(
        {
          values$query_geno_react <- query_genotypes(
            db_path = input$db_path,
            variant_ids = query_by_alf_result()$variant_id,
            meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
          )

          show_toast_success(text = "Putative Causal Variants Found!")
        },
        error = function(e) {
          show_alert(
            title = "No Putative Causal Variants Found!",
            text = "Confirm Impact Level and Allele Frequency Threshold Using Annotation Summary",
            type = "error",
            showCloseButton = TRUE,
            timer = 5000
          )
        }
      )

      # Display result for successful pcvs
      output$pcvs_kasp_marker_design_result <- renderUI({
        genotype_results_ui(ns)
      })
    })


    # Result genotype matrix
    output$genotype_results_tbl <- DT::renderDT({
      req(values$query_geno_react)
      render_dt(values$query_geno_react)
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

      showModal(
        modalDialog(
          title = div(
            style = "text-align: center;",
            "Design KASP Marker From Putative Causal Variants"
          ),
          size = "xl",
          footer = modalButton("Close"),
          fluidRow(
            column(
              width = 4,
              bslib::card(
                bslib::card_header("Analysis Parameters",
                  class = "bg-primary text-white"
                ),
                bslib::card_body(
                  fileInput(
                    ns("modal_genome_file"),
                    label = "Genome Reference File",
                    accept = c(".fa", ".fasta", ".gz")
                  ),
                  selectizeInput(
                    ns("modal_marker_ID"),
                    label = "Marker ID",
                    choices = values$query_geno_react$variant_id,
                    options = list(placeholder = "Select a variant...")
                  ),
                  textInput(
                    ns("modal_reg_name"),
                    label = "Region Name",
                    placeholder = "e.g., drought resistance locus"
                  ),
                  numericInput(
                    ns("modal_maf"),
                    label = "Minor Allele Frequency (MAF)",
                    value = 0.05, min = 0, max = 1, step = 0.01
                  ),
                  bslib::input_switch(
                    ns("modal_draw_plot"),
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
                    class = "btn-primary"
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
                  bslib::card(
                    bslib::card_footer(
                      fluidRow(
                        column(width = 3, selectInput(
                          inputId = ns("exten"),
                          label = "Download file as?",
                          choices = c(".csv", ".xlsx"),
                          selected = ".csv",
                          multiple = FALSE
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
                  )
                ),
                uiOutput(ns("plot_container"))
              )
            )
          )
        )
      )
    })

    # Server Side of kasp marker design
    kasp_design_result <- reactiveVal(NULL) # store kasp result

    observeEvent(input$modal_run_but, {
      req(
        values$query_geno_react, input$modal_genome_file$datapath,
        input$modal_maf, input$modal_marker_ID
      )

      shinyWidgets::show_toast(
        title = "Designing Marker...",
        text = "Please Wait...",
        timer = 8000,
        type = "info"
      )

      tryCatch(
        {
          # Get column names
          gt_cols <- colnames(values$query_geno_react)

          # Find column indices
          id_col <- gt_cols[grep("id", gt_cols, ignore.case = TRUE)[1]]
          chrom_col <- gt_cols[grep("chro", gt_cols, ignore.case = TRUE)[1]]
          pos_col <- gt_cols[grep("pos", gt_cols, ignore.case = TRUE)[1]]
          ref_col <- gt_cols[grep("ref", gt_cols, ignore.case = TRUE)[1]]
          alt_col <- gt_cols[grep("alt", gt_cols, ignore.case = TRUE)[1]]

          # Hard coded genotype start.
          geno_start <- 7

          # Get unique chromosomes
          unique_chr <- unique(values$query_geno_react[[chrom_col]])

          # Run KASP design
          result_data <- kasp_marker_design_alt(
            vcf_file = NULL,
            gt_df = values$query_geno_react,
            variant_id_col = id_col,
            chrom_col = chrom_col,
            pos_col = pos_col,
            ref_al_col = ref_col,
            alt_al_col = alt_col,
            geno_start = geno_start,
            marker_ID = input$modal_marker_ID,
            chr = unique_chr,
            genome_file = input$modal_genome_file$datapath,
            plot_file = tempdir(),
            region_name = input$modal_reg_name,
            maf = input$modal_maf,
            plot_draw = TRUE
          )

          kasp_design_result(result_data)

          # Show success message
          shinyWidgets::show_toast(
            title = "Success",
            text = "KASP marker designed successfully",
            type = "success",
            timer = 5000
          )
        },
        error = function(e) {
          # Show error message
          shinyWidgets::show_toast(
            title = "Error",
            text = paste("Failed to design KASP marker:", e$message),
            type = "error",
            timer = 8000
          )
          kasp_design_result(NULL)
        }
      )
    })

    output$kasp_table <- DT::renderDT({
      req(kasp_design_result()$marker_data)
      render_dt(data = kasp_design_result()$marker_data)
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
          write.csv(kasp_design_result()$marker_data, file, row.names = FALSE)
        } else if (input$exten == ".xlsx") {
          openxlsx::write.xlsx(kasp_design_result()$marker_data, file)
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

    # Render alignment Plot
    output$plot <- renderPlot({
      req(kasp_design_result())

      print(kasp_design_result()$plot) # show plot
    })

    # Download Plot
    output$download_plot <- downloadHandler(
      filename = function() {
        # Clean proposed user file name and append underscores
        clean_marker <- gsub("[^[:alnum:]_-]", "_", input$modal_marker_ID)
        paste0("alignment_", clean_marker, ".pdf")
      },
      content = function(file) {
        ggplot2::ggsave(
          filename = file,
          plot = kasp_design_result()$plot,
          device = "pdf",
          width = 24,
          height = 9,
          units = "in"
        )
      }
    )
  })
}

## To be copied in the UI
# mod_variant_discovery_ui("variant_discovery_1")

## To be copied in the server
# mod_variant_discovery_server("variant_discovery_1")
