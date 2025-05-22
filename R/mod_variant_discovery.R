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
#'
#'
mod_variant_discovery_ui <- function(id) {
  ns <- NS(id)
  tagList(
    bslib::navset_card_underline(
      shinyjs::useShinyjs(),
      sidebar = bslib::sidebar(
        width = 350,
        title = "Database Connection",
        status = "primary",
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
              #   h4("Database Overview"),
              verbatimTextOutput(ns("db_info")) # text output for database.
            )
          )
        )
      ),

      # Variant Discovery Tab
      bslib::nav_panel(
        title = tags$b("Info"), icon = icon("info-circle"),
        conditionalPanel(
          condition = paste0("output['", ns("is_connected"), "'] == true"),
          # Variant impact table
          bslib::layout_column_wrap(
            width = "100%",
            bslib::card(
              bslib::card_header(
                class = "bg-primary text-white",
                "Variant Impact Summary"
              ),
              bslib::card_body(
                DT::dataTableOutput(ns("table_impact_id"))
              )
            )
          ),
          bslib::layout_column_wrap(
            width = "100%",
            bslib::card(
              bslib::card_header(
                class = "bg-primary text-white",
                "Variant Statistics"
              ),
              bslib::card_body(
                DT::dataTableOutput(ns("table_var_stats_id"))
              )
            )
          ),
          bslib::layout_column_wrap(
            width = 1 / 2,
            bslib::card(
              bslib::card_header(
                class = "bg-primary text-white",
                "Summarised SQLite Tables"
              ),
              bslib::card_body(
                DT::dataTableOutput(outputId = ns("sum_sqlite_id"))
              )
            ),
            bslib::card(
              bslib::card_header(
                class = "bg-primary text-white",
                "Variant Type Count"
              ),
              bslib::card_body(
                DT::dataTableOutput(outputId = ns("count_variant_typ_id"))
              )
            )
          ),
          bslib::layout_column_wrap(
            width = "100%",
            bslib::card(
              bslib::card_header(
                class = "bg-primary text-white",
                "Inspect SQlite Table Schema"
              ),
              bslib::card_body(
                selectInput(
                  inputId = ns("table_name_lst"),
                  label = "Table to Query",
                  choices = c("variants", "annotations", "genotypes"),
                  selected = "genotypes"
                ),
                DT::dataTableOutput(outputId = ns("results_lst"))
              )
            )
          )
        )
      ),

      # Variant Analysis Tab
      bslib::nav_panel(
        title = tags$b("Query Actions"), icon = icon("bolt"),
        # icon = icon("chart-bar"),
        conditionalPanel(
          condition = paste0("output['", ns("is_connected"), "'] == false"),
          bslib::card(
            bslib::card_header("Variant Analysis"),
            bslib::card_body(
              "Please connect to a database to use this feature."
            )
          )
        ),
        conditionalPanel(
          condition = paste0("output['", ns("is_connected"), "'] == true"),

          # Show fixed card for Query database.
          bslib::layout_column_wrap(
            width = 1 / 2,
            bslib::card(
              bslib::card_header(tags$strong("Query Database")),
              # Modal show if user does not have the genocordinates.
              helpText(h5("Don't have genotype co-ordinates?")),
              div(
                style = "display: flex; justify-content: center;",
                actionButton(
                  inputId = ns("get_cord"), label = "Get Genotype Co-ordinates",
                  width = "70%"
                )
              ),
              selectInput(
                inputId = ns("table_name"),
                label = "Table to Query",
                choices = c("variants", "annotations", "genotypes"),
                selected = "annotations"
              ),
              textInput(
                inputId = ns("chrom"),
                label = "Chromosome",
                value = NULL
              ),
              numericInput(
                inputId = ns("start"),
                label = "Start Position",
                value = NULL
              ),
              numericInput(
                inputId = ns("end"),
                label = "End Position",
                value = NULL
              ),
              uiOutput(outputId = ns("gene_name_id")),
              bslib::card_footer(
                actionButton(
                  inputId = ns("query_db_btn"),
                  label = "Query Database",
                  class = "btn-primary w-100",
                  icon = icon("database")
                )
              ),
              style = "max-width: 350px;"
            ),
            # Section for rendering based on users input.
            bslib::card(
              bslib::card_header(selectInput(
                inputId = ns("choose_query"),
                label = span("Choose a query Action", style = "margin-right: 10px;"),
                multiple = FALSE,
                choices = c(
                  "None",
                  "Query Annotation Summary",
                  "Query by Impact",
                  "Query by Allele Frequency",
                  "Query by Genotype"
                )
              )),
              uiOutput(ns("result_display")),
              style = "max-width: 400px;"
            )
          ),
          # uiOutput(ns("shiny_mod")), # show modal UI

          bslib::navset_card_tab(
            id = ns("nav_id"), selected = "Main Database Query Result",
            bslib::nav_panel(
              title = "Main Database Query Result",
              uiOutput(outputId = ns("query_db_display"))
            ),
            bslib::nav_panel(
              title = "Other Query Action Result",
              uiOutput(outputId = ns("other_db_result")) # display query by annotation summary
            ),
            # bslib::nav_panel(
            #   title = "Query by Impact",
            #   uiOutput(outputId = ns("query_by_impact_result")) # display query by impact. and user can filter
            # ),
            # bslib::nav_panel(
            #   title = "Filter by allele frequency",
            #   uiOutput(outputId = ns("filter_by_allele_freq_result")) # Filter by allele frequency
            # ),
            # bslib::nav_panel(
            #   title = "Query by Genotype",
            #   uiOutput(outputId = ns("query_by_genotype_result")) # Query by genotype
            # )
          )
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
      connected = FALSE, # Connection status
      tables = NULL, # List of database tables
      status_message = "Not connected", # Status message
      selected_table = NULL, # Currently selected table
      variant_impact = NULL, # Variant impact analysis results
      db_path = NULL, # Current database path
      sqlite_summary = NULL, # store sqlite summary
      variant_count = NULL, # store variant count
      variant_stats = NULL # store variant statistics
    )


    # Display status (Connected or not connected)
    output$status_badge <- renderUI({
      if (rv$connected) {
        span(
          class = "status-connected", style = "margin-left: 10px;",
          icon("check-circle"), " Connected"
        )
      } else {
        span(
          class = "status-disconnected", style = "margin-left: 10px;",
          icon("times-circle"), " Disconnected"
        )
      }
    })

    #--------------------------------------------
    # DATABASE CONNECTION MANAGEMENT
    #--------------------------------------------
    # Connect to database
    observeEvent(input$connect_btn, {
      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Connecting to Database... Please wait."
      )
      tryCatch({
        req(input$db_path)


        # Validate the path exists
        if (!file.exists(input$db_path)) {
          showNotification("Error: Database file not found", type = "error")
          return()
        }

        # Close existing connection if any
        if (!is.null(conn)) {
          dbDisconnect(conn)
          conn <- NULL
        }

        # Connect using the CLEAN path
        conn <<- dbConnect(RSQLite::SQLite(), input$db_path)

        # Get tables using clean path
        tables <- list_sqlite_tables(input$db_path)

        # Update reactive values with clean path
        rv$variant_impact <- variant_impact_summary(db_path = input$db_path)
        rv$sqlite_summary <- summarize_sqlite_tables(db_path = input$db_path)
        rv$variant_count <- count_variant_types(db_path = input$db_path)
        rv$variant_stats <- variant_stats(db_path = input$db_path)

        # Update state
        rv$connected <- TRUE
        rv$tables <- tables
        rv$status_message <- paste("Connected with", length(tables), "tables")
        rv$db_path <- input$db_path # Store normalized path

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


    #-----------
    # If user disconnects the database.
    #-----------
    # Disconnect from database
    observeEvent(input$disconnect_btn, {
      if (!is.null(conn)) {
        dbDisconnect(conn)
        conn <<- NULL

        # Update state
        rv$connected <- FALSE
        rv$tables <- NULL
        rv$status_message <- "Disconnected"
        rv$variant_impact <- NULL
        rv$selected_table <- NULL

        # Show toast
        shinyWidgets::show_toast(title = 'Disconnected from Database',
                                 type = 'warning',
                                 timer = 5000)
      }
    })

    # Clean up connection when app closes
    onSessionEnded(function() {
      if (!is.null(conn)) {
        dbDisconnect(conn)
      }
    })

    # Connection status for conditional panels
    output$is_connected <- reactive({
      rv$connected
    })
    outputOptions(output, "is_connected", suspendWhenHidden = FALSE)

    # Database information output
    output$db_info <- renderPrint({
      req(rv$connected)

      # Get basic database info
      cat("Path: ", rv$db_path, "\n")
      cat("Tables: ", length(rv$tables), "\n")
      cat("Available tables:\n")
      cat(paste("- ", rv$tables, collapse = "\n"))
    })

    #---------------------
    # Info section once database is loaded
    #------------------
    # Output variant impact
    output$table_impact_id <- renderDT({
      req(rv$variant_impact)
      datatable(rv$variant_impact)
    })

    # Output variant statistics
    output$table_var_stats_id <- renderDT({
      req(rv$variant_stats)
      datatable(rv$variant_stats)
    })

    # Output sqlite table summary
    output$sum_sqlite_id <- renderDT({
      req(rv$sqlite_summary)
      datatable(rv$sqlite_summary)
    })

    # Output variant count
    output$count_variant_typ_id <- renderDT({
      req(rv$variant_count)
      datatable(rv$variant_count)
    })

    # list tables columns
    lst_tbl_column <- reactive({
      req(input$table_name_lst, input$db_path)
      list_table_columns(
        db_path = input$db_path,
        table_name = input$table_name_lst
      )
    })

    # Output.
    output$results_lst <- renderDT({
      req(lst_tbl_column())
      datatable(lst_tbl_column())
    })

    #-------------
    # Query section, UI rendering.
    #-------------

    # If user does not have genotype cordinates to input manually.
    observeEvent(input$get_cord, {
      showModal(
        modalDialog(
          easyClose = FALSE,
          card(
            card_header("Gene Parameters"),
            textInput(
              inputId = ns("gene_name"), label = "Gene Name (Sobic ID)",
              value = "Sobic.005G213600",
              placeholder = "Enter Sobic ID (e.g., Sobic.005G213600)"
            ),
            radioButtons(
              inputId = ns("input_method"), label = "GFF File Source",
              choices = c("URL" = "url", "Upload File" = "file"),
              selected = "url"
            ),
            conditionalPanel(
              condition = paste0("input['", ns("input_method"), "'] === 'url'"),
              textInput(
                inputId = ns("gff_url"), label = "GFF File URL",
                value = "https://raw.githubusercontent.com/awkena/panGB/main/Sbicolor_730_v5.1.gene.gff3.gz",
                placeholder = "Enter URL to GFF file"
              )
            ),
            conditionalPanel(
              condition = paste0("input['", ns("input_method"), "'] === 'file'"),
              fileInput(
                inputId = ns("gff_file"), label = "Upload GFF File",
                accept = c(".gff3", ".gff", ".gff3.gz", ".gff.gz")
              )
            ),
            card_footer(
              actionButton(
                inputId = ns("submit"), label = "Get Coordinates",
                class = "btn-primary w-100",
                icon = icon("search")
              )
            )
            # style = "max-width: 350px;"
          )
        )
      )
    })


    # Close modal window.
    observeEvent(input$submit, {
      removeModal()

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Getting gene coordinates... Please wait."
      )
      # Get the gene name
      gene_name <- input$gene_name

      # Determine GFF path based on input method
      if (input$input_method == "url") {
        gff_path <- input$gff_url
      } else {
        gff_path <- input$gff_file$datapath
      }

      # Attempt to get gene coordinates
      tryCatch({
        result <- gene_coord_gff(gene_name, gff_path)
        values$result <- result

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

    # Populate all the genotype sections.
    observeEvent(values$result, {
      # All chromosome section
      updateTextInput(session, inputId = "chrom", value = values$result$chrom)
      updateNumericInput(session, inputId = "start", value = values$result$start)
      updateNumericInput(session, inputId = "end", value = values$result$end)
    })

    # Update ui for database.
    observeEvent(input$query_db_btn, {
      updateTabsetPanel(session, inputId = "nav_id", selected = "Main Database Query Result")

      output$query_db_display <- renderUI({
        tagList(
          # Query database table result...
          card(
            card_header("Query Results"),
            div(
              style = "overflow-x: auto;",
              DT::DTOutput(outputId = ns("query_results"))
            ),
            bslib::card_footer(
              # File naming input
              textInput(
                inputId = ns("File_name"),
                label = "Enter File Name",
                value = "Chrom 05",
                width = "30%"
              ),

              # Action buttons container
              div(
                style = "display: flex; justify-content: space-between; align-items: center; gap: 10px; width: 100%;",

                # Download Excel button
                downloadButton(
                  outputId = ns("download_excel"),
                  label = "Export to Excel",
                  class = "btn btn-primary",
                  style = "flex-grow: 1;"
                ),

                # Push to KASP button
                actionButton(
                  inputId = ns("push"),
                  label = "Push to KASP Marker Design Tab",
                  class = "btn-success",
                  style = "flex-grow: 1; white-space: normal;"
                )
              )
            ),
            class = "mt-3"
          )
          # conditionalPanel(
          #   condition = paste0("input['", ns("table_name"), "'] === 'genotypes'"),
          #   input_switch(id = ns("filter_switch"), label = "Filter table", value = FALSE),
          #   conditionalPanel(
          #     condition = paste0("input['", ns("filter_switch"), "'] === true"),
          #     card(
          #       card_header("Filter by Allele Frequency"),
          #       fluidRow(
          #         column(
          #           3,
          #           textInput(inputId = ns("variant_id_col"), label = "Variant ID Column", value = "variant_id")
          #         ),
          #         column(
          #           3,
          #           textInput(inputId = ns("chrom_col"), label = "Chromosome Column", value = "chrom")
          #         ),
          #         column(
          #           3,
          #           textInput(inputId = ns("pos_col"), label = "Position Column", value = "pos")
          #         ),
          #         column(
          #           3,
          #           sliderInput(
          #             inputId = ns("filter_range"), label = "Allele Frequency Range",
          #             min = 0, max = 1, value = c(0, 1), step = 0.01
          #           )
          #         )
          #       ), DT::DTOutput(outputId = ns("filter_result_tbl")),
          #       class = "mt-3"
          #     ),
          #     input_switch(id = ns("calc_af"), label = "Compute Allele Frequency", value = FALSE)
          #   ),
          #   conditionalPanel(
          #     condition = paste0("input['", ns("calc_af"), "'] === true"),
          #     card(
          #       card_header("Computed Allele Frequency"),
          #       fluidRow(
          #         column(
          #           3,
          #           textInput(inputId = ns("variant_id_col_af"), label = "Variant ID Column", value = "variant_id")
          #         ),
          #         column(
          #           3,
          #           textInput(inputId = ns("chrom_col_af"), label = "Chromosome Column", value = "chrom")
          #         ),
          #         column(
          #           3,
          #           textInput(inputId = ns("pos_col_af"), label = "Position Column", value = "pos")
          #         )
          #       ), DT::DTOutput(outputId = ns("calc_af_result_tbl"))
          #     ),
          #     class = "mt-3"
          #   )
          # )
        )
      })
    })

    #----------------
    # Other database queries ui
    #----------------

    # Observer for select inputs
    observeEvent(input$choose_query, {
      req(input$choose_query)

      if (input$choose_query == "Query by Impact") {
        output$result_display <- renderUI({
          tagList(
            card(
              card_header("Impact Parameters"),
              selectInput(
                inputId = ns("impact_level"), label = "Impact Level",
                choices = c("HIGH", "MODERATE", "LOW", "MODIFIER")
              ),
              textInput(
                inputId = ns("impact_chrom"),
                label = "Chromosome",
                value = if (!is.null(values$result)) values$result$chrom else NULL
              ),
              numericInput(
                inputId = ns("impact_start"),
                "Start Position",
                value = if (!is.null(values$result)) values$result$start else NULL
              ),
              numericInput(
                inputId = ns("impact_end"),
                "End Position",
                value = if (!is.null(values$result)) values$result$end else NULL
              ),
              card_footer(
                actionButton(
                  inputId = ns("query_impact_btn"), label = "Query by Impact",
                  class = "btn-primary w-100",
                  icon = icon("filter")
                )
              ),
              style = "max-width: 350px;"
            )
          )
        })

        output$other_db_result <- renderUI({
          card(
            card_header("Variants by Impact"),
            div(
              style = "overflow-x: auto;",
              DT::DTOutput(outputId = ns("impact_results"))
            ),
            class = "mt-3"
          )
        })
      } else if (input$choose_query == "Query by Allele Frequency") {
        output$result_display <- renderUI({
          tagList(
            card(
              card_header("Allele Frequency Parameters"),
              sliderInput(
                inputId = ns("af_range"), label = "Allele Frequency Range",
                min = 0, max = 1, value = c(0, 1), step = 0.01
              ),
              textInput(
                inputId = ns("af_chrom"),
                label = "Chromosome ",
                value = if (!is.null(values$result)) values$result$chrom else NULL
              ),
              numericInput(
                inputId = ns("af_start"), label = "Start Position",
                value = if (!is.null(values$result)) values$result$start else NULL
              ),
              numericInput(
                inputId = ns("af_end"), label = "End Position",
                value = if (!is.null(values$result)) values$result$end else NULL
              ),
              card_footer(
                actionButton(
                  inputId = ns("query_af_btn"), label = "Filter by AF",
                  class = "btn-primary w-100",
                  icon = icon("filter")
                )
              ),
              style = "max-width: 350px;"
            )
          )
        })

        output$other_db_result <- renderUI({
          card(
            card_header("Variants Filtered by AF"),
            div(
              style = "overflow-x: auto;",
              DT::DTOutput(outputId = ns("af_results"))
            ),
            class = "mt-3"
          )
        })
      } else if (input$choose_query == "Query by Genotype") {
        output$result_display <- renderUI({
          tagList(
            card(
              card_header("Variant Selection"),
              textAreaInput(
                inputId = ns("variant_ids"), label = "Variant IDs (one per line)",
                height = "150px"
              ),
              card_footer(
                actionButton(
                  inputId = ns("query_genotypes_btn"), label = "Get Genotypes",
                  class = "btn-primary w-100",
                  icon = icon("dna")
                )
              ),
              style = "max-width: 350px;"
            )
          )
        })

        output$other_db_result <- renderUI({
          card(
            card_header("Genotype Data"),
            div(
              style = "overflow-x: auto;",
              DT::DTOutput(outputId = ns("genotype_results_tbl"))
            ),
            class = "mt-3",
            bslib::card_footer(
              # File naming input
              textInput(
                inputId = ns("File_name"),
                label = "Enter File Name",
                value = "Chrom 05",
                width = "30%"
              ),
              # Action buttons container
              div(
                style = "display: flex; justify-content: space-between; align-items: center; gap: 10px; width: 100%;",

                # Download Excel button
                downloadButton(
                  outputId = ns("download_excel"),
                  label = "Export to Excel",
                  class = "btn btn-primary",
                  style = "flex-grow: 1;"
                ),

                # Push to KASP button
                actionButton(
                  inputId = ns("push"),
                  label = "Push to KASP Marker Design Tab",
                  class = "btn-success",
                  style = "flex-grow: 1; white-space: normal;"
                )
              )
            )
          )
        })
      } else if (input$choose_query == "Query Annotation Summary") {
        output$result_display <- renderUI({
          tagList(
            card(
              card_header("Region Parameters"),
              textInput(
                inputId = ns("ann_chrom"), label = "Chromosome",
                value = if (!is.null(values$result)) values$result$chrom else NULL
              ),
              selectInput(
                inputId = ns("table_name_a"),
                label = "Table name for Annotation",
                choices = c("variants", "annotations", "genotypes"),
                selected = "annotations"
              ),
              selectInput(
                inputId = ns("table_name_v"),
                label = "Table name for Variants",
                choices = c("variants", "annotations", "genotypes"),
                selected = "variants"
              ),
              numericInput(
                inputId = ns("ann_start"), label = "Start Position",
                value = if (!is.null(values$result)) values$result$start else NULL
              ),
              numericInput(
                inputId = ns("ann_end"), label = "End Position",
                value = if (!is.null(values$result)) values$result$end else NULL
              ),
              card_footer(
                actionButton(
                  inputId = ns("get_ann_summary_btn"), "Get Summary",
                  class = "btn-primary w-100",
                  icon = icon("chart-pie")
                )
              ),
              style = "max-width: 350px;"
            )
          )
        })

        output$other_db_result <- renderUI({
          card(
            card_header("Annotation Statistics"),
            card(
              card_header("Annotation Summary", class = "bg-primary text-light"),
              DT::DTOutput(outputId = ns("ann_summary_tbl"))
            ),
            card(
              card_header("Impact Summary", class = "bg-primary text-light"),
              DT::DTOutput(outputId = ns("impact_summary_tbl"))
            ),
            card(
              card_header("Variant Type Totals", class = "bg-primary text-light"),
              DT::DTOutput(outputId = ns("variant_totals_tbl"))
            ),
            class = "mt-3",
          )
        })
      } else {
        output$result_display <- renderUI({
          NULL
        })
      }
    })


    # Empty reactive values for use
    values <- reactiveValues(
      status = "",
      result = NULL,
      query_impact_react = NULL,
      query_af_react = NULL,
      query_db_val = NULL,
      query_ann_react = NULL,
      query_geno_react = NULL
    )

    # Query database--
    # Render UI for gene name if user selects annotation
    observeEvent(input$table_name, {
      req(input$table_name)
      if (input$table_name == "annotations") {
        output$gene_name_id <- renderUI({
          textInput(
            inputId = ns("query_gene_name"),
            label = "Gene Name",
            value = NULL
          )
        })
      } else {
        output$gene_name_id <- renderUI(NULL)
      }
    })

    # Observe when user clicks to query database.
    observeEvent(input$query_db_btn, {
      req(input$start, input$end, input$chrom, input$db_path, input$table_name)
      shinybusy::show_modal_spinner(
        spin = "fading-circle", # Choose a spinner style (e.g., "atom", "circle")
        color = "#0dc5c1", # Spinner color
        text = "Querying Database... Please wait."
      )
      tryCatch({
        values$query_db_val <- query_db(
          db_path = input$db_path,
          table_name = input$table_name,
          chrom = input$chrom,
          start = input$start,
          end = input$end,
          gene_name = if (input$query_gene_name == "") NULL else input$query_gene_name
        )

        show_alert(
          title = "Success!",
          # text = "Found Gene Cordinates",
          type = "success",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, error = function(e) {
        show_alert(
          title = "Failed!",
          text = "Couln't Query Database",
          type = "danger",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })


    #----------------------------
    # Rendering of results section
    #---------------------------

    output$query_results <- renderDT({
      req(values$query_db_val)
      datatable(values$query_db_val)
    })

    # Download result as excel sheet.
    output$download_excel <- downloadHandler(
      filename = function() {
        paste("dbquery", input$File_name, input$table_name, ".xlsx", sep = "_")
      },
      content = function(file) {
        writexl::write_xlsx(values$query_db_val, path = file)
      }
    )


    # filter section---
    filter_by_af_react <- reactive({
      req(
        values$query_db_val, input$variant_id_col,
        input$chrom_col, input$pos_col, input$filter_range
      )
      filter_by_af(
        gt = values$query_db_val,
        variant_id_col = input$variant_id_col,
        chrom_col = input$chrom_col,
        pos_col = input$pos_col,
        min_af = input$filter_range[1],
        max_af = input$filter_range[2]
      )
    })

    output$filter_result_tbl <- renderDT({
      req(filter_by_af_react())
      datatable(filter_by_af_react())
    })

    # # Calculate allele frequency---
    # calc_af_result <- reactive({
    #   req(values$query_db_val, input$variant_id_col_af, input$chrom_col_af, input$pos_col_af)
    #   panGenomeBreedr::calc_af(
    #     gt = values$query_db_val,
    #     variant_id_col = input$variant_id_col_af,
    #     chrom_col = input$chrom_col_af,
    #     pos_col = input$pos_col_af
    #   )
    # })
    #
    # output$calc_af_result_tbl <- renderDT({
    #   req(calc_af_result())
    #   datatable(calc_af_result())
    # })

    # Query by impact---
    observeEvent(input$query_impact_btn, {
      updateTabsetPanel(session, inputId = "nav_id", selected = "Other Query Action Result")
      req(
        input$db_path, input$impact_level,
        input$impact_chrom, input$impact_start, input$impact_end
      )

      # Stores result for query by impact
      values$query_impact_react <- query_by_impact(
        db_path = input$db_path,
        impact_level = input$impact_level,
        chrom = input$impact_chrom,
        start = input$impact_start,
        end = input$impact_end
      )
      # Show toast
      shinyWidgets::show_toast(
        title = "Success",
        text = "Variants Filtered By Impact Level",
        timer = 2000,
        position = "bottom-end",
        type = "success"
      )
    })
    output$impact_results <- renderDT({
      req(values$query_impact_react)
      datatable(values$query_impact_react)
    })


    # Get genotypes selected by the individual based on impact.
    hold_genotypes_impact <- reactive({
      req(values$query_impact_react, input$db_path) # need

      query_genotypes(
        db_path = input$db_path,
        variant_ids = values$query_impact_react$variant_id,
        meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
      )
    })



    # Query allele frequency---
    observeEvent(input$query_af_btn, {
      req(
        input$af_range,
        input$af_chrom, input$af_start, input$af_end,
        hold_genotypes_impact()
      )
      updateTabsetPanel(session, inputId = "nav_id", selected = "Other Query Action Result")

      values$query_af_react <- filter_by_af(
        gt = hold_genotypes_impact(),
        min_af = input$af_range[1],
        max_af = input$af_range[2]
      )
      # Show toast
      shinyWidgets::show_toast(
        title = "Success",
        text = "Variants Filtered By Allele Frequency",
        timer = 2000,
        position = "bottom-end",
        type = "success"
      )
    })

    output$af_results <- renderDT({
      req(values$query_af_react)
      datatable(values$query_af_react)
    })



    # Query by genotypes---
    observeEvent(input$query_genotypes_btn, {
      updateTabsetPanel(session, inputId = "nav_id", selected = "Other Query Action Result")
      req(values$query_impact_react, input$db_path) # need

      values$query_geno_react <- query_genotypes(
        db_path = input$db_path,
        variant_ids = values$query_af_react$variant_id,
        meta_data = c("chrom", "pos", "ref", "alt", "variant_type")
      )
      # Show toast
      shinyWidgets::show_toast(
        title = "Success",
        text = "Variant Filtered By Genotype Calls",
        timer = 2000,
        position = "bottom-end",
        type = "success"
      )
    })

    output$genotype_results_tbl <- renderDT({
      req(values$query_geno_react)
      datatable(values$query_geno_react)
    })

    # Query annotation summary---
    observeEvent(input$get_ann_summary_btn, {
      updateTabsetPanel(session, inputId = "nav_id", selected = "Other Query Action Result")
      req(
        input$db_path, input$ann_end, input$table_name_v,
        input$ann_chrom, input$ann_start, input$table_name_a
      )

      values$query_ann_react <- query_ann_summary(
        db_path = input$db_path,
        variants_table = input$table_name_v,
        annotations_table = input$table_name_a,
        chrom = input$ann_chrom,
        start = input$ann_start,
        end = input$ann_end
      )
      # Show toast
      shinyWidgets::show_toast(
        title = "Success",
        text = "Variant Filtered By Annotation Summary",
        timer = 2000,
        position = "bottom-end",
        type = "success"
      )
    })

    output$ann_summary_tbl <- renderDT({
      req(values$query_ann_react)
      datatable(values$query_ann_react$annotation_summary)
    })

    output$impact_summary_tbl <- renderDT({
      req(values$query_ann_react)
      datatable(values$query_ann_react$impact_summary)
    })

    output$variant_totals_tbl <- renderDT({
      req(values$query_ann_react)
      datatable(values$query_ann_react$variant_type_totals)
    })

    # Reactive value to store data to return
    return_value <- reactiveVal(NULL)

    # Observe the push button
    observeEvent(input$push, {
      req(values$query_geno_react)

      df <- as.data.frame(values$query_geno_react)
      return_value(df)

      shinyWidgets::show_toast(
        title = "Success",
        text = "Proceed to KASP Marker Design",
        type = "success",
        timer = 5000
      )
    })

    # Return to be used in kasp marker
    return(reactive(return_value()))
  })
}

## To be copied in the UI
# mod_variant_discovery_ui("variant_discovery_1")

## To be copied in the server
# mod_variant_discovery_server("variant_discovery_1")
