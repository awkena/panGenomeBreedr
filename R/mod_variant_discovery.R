#' Variant Discovery UI Function
#'
#' @description A shiny Module for variant discovery that allows database connections
#'   and various query operations for genomic variants.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @return A UI definition
#'
#' @noRd
#'
#' @importFrom shiny NS tagList conditionalPanel textInput actionButton
#' @importFrom shiny uiOutput selectInput numericInput verbatimTextOutput
#' @importFrom shiny icon div h4 span helpText tags
#' @importFrom bslib navset_card_underline sidebar card card_header
#' @importFrom bslib card_body card_footer nav_panel navset_card_tab
#' @importFrom bslib layout_column_wrap input_switch
#' @importFrom shinyjs useShinyjs
#' @importFrom DT dataTableOutput
#' @importFrom dplyr %>%
mod_variant_discovery_ui <- function(id) {
  ns <- NS(id)

  # Reusable card styling
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
          bslib::card_body(DT::dataTableOutput(ns("table_impact_id")))
        ),
        bslib::card(
          primary_card_header("Variant Statistics"),
          bslib::card_body(DT::dataTableOutput(ns("table_var_stats_id")))
        )
      ),
      bslib::layout_column_wrap(
        width = 1 / 2,
        bslib::card(
          primary_card_header("Summarised SQLite Tables"),
          bslib::card_body(DT::dataTableOutput(ns("sum_sqlite_id")))
        ),
        bslib::card(
          primary_card_header("Variant Type Count"),
          bslib::card_body(DT::dataTableOutput(ns("count_variant_typ_id")))
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
            DT::dataTableOutput(ns("results_lst"))
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
          bslib::card_header("Variant Analysis"),
          bslib::card_body("Please connect to a database to use this feature.")
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
        value = c(0, 1),
        step = 0.01,
        width = "100%"
      ),
      bslib::card_footer(
        actionButton(
          inputId = ns("get_pcv_btn"), label = "Exract PCVs",
          class = "btn-primary",
          width = "100%",
          icon = icon("filter")
        )
      )
    )
  }

  # Main UI Structure
  tagList(
    bslib::navset_card_underline(
      shinyjs::useShinyjs(),
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
        info_tab(ns)
      ) # Display info
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
#' @importFrom DT dataTableOutput DTOutput
#'
#' @noRd
mod_variant_discovery_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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

    #--------------------------------------------
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
          showNotification("Error: Database file not found", type = "error")
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

    #---------------------
    # Info section once database is loaded
    #------------------

    # Users decide to get info about table then it renders.
    observeEvent(input$get_info, {
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Retrieving database information… This may take some time"
      )

      tryCatch({
        req(input$db_path)

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

    # List table columns
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

    #-------------
    # Query section, UI rendering.
    #-------------
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

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Getting gene coordinates... Please wait."
      )

      tryCatch({
        gff_path <- if (input$input_method == "url") input$gff_url else input$gff_file$datapath
        values$result <- gene_coord_gff(input$gene_name, gff_path)

        show_alert(
          title = "Success!",
          text = "Found Gene Coordinates",
          type = "success",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, error = function(e) {
        show_alert(
          title = "Failed!",
          text = "Couldn't find Gene Coordinates",
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })



    # If user selects to manually input.
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
          text = "Couln't Query Database,Set Genotype Co-ordinate first!",
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })

    ######################### Annotation summary.
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



    # Download handler
    output$download_excel <- downloadHandler(
      filename = function() {
        paste("dbquery", input$File_name, input$table_name, ".xlsx", sep = "_")
      },
      content = function(file) {
        writexl::write_xlsx(values$query_db_val, path = file)
      }
    )

    #----------------
    # Other database queries
    #----------------
    # This is the only UI I will have to use.
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
            style = "display: flex; justify-content: space-between; align-items: center; gap: 10px; width: 100%;",
            downloadButton(
              outputId = ns("download_excel"),
              label = "Export to Excel",
              class = "btn btn-primary",
              style = "flex-grow: 1;"
            ),
            actionButton(
              inputId = ns("push"),
              label = "Push to KASP Marker Design Tab",
              class = "btn-success",
              style = "flex-grow: 1; white-space: normal;"
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
      tryCatch(
        {
          req(query_by_alf_result(), input$db_path)

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

    # Return value for KASP marker
    return_value <- reactiveVal(NULL)

    observeEvent(input$push, {
      req(values$query_geno_react)
      return_value(as.data.frame(values$query_geno_react))

      shinyWidgets::show_toast(
        title = "Success",
        text = "Proceed to KASP Marker Design",
        type = "success",
        timer = 5000
      )
    })

    return(reactive(return_value()))
  })
}

## To be copied in the UI
# mod_variant_discovery_ui("variant_discovery_1")

## To be copied in the server
# mod_variant_discovery_server("variant_discovery_1")
