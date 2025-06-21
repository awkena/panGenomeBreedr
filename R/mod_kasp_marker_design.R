#' KASP Marker Design UI Function
#'
#' @description A Shiny module for KASP marker design
#'
#' @param id Internal parameter for {shiny}.
#'
#' @return A UI for KASP Marker Design
#'
#' @noRd
#'
#' @importFrom shiny NS tagList radioButtons
#' @importFrom bslib layout_sidebar sidebar card card_header card_footer input_switch accordion accordion_panel
#'
mod_kasp_marker_design_ui <- function(id) {
  ns <- NS(id)
  tagList(
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = "KASP Marker Design Parameters",
        width = 400,

        # File Uploads Card
        bslib::card(
          bslib::card_header(tags$b("File Uploads")),
          fileInput(ns("genome_file"),
            label = "Genome File (.fa, .fasta, .gz)",
            accept = c(".fa", ".fasta", ".gz")
          ),
          radioButtons(
            inputId = ns("upload_choice"),
            label = "Variant Data Type",
            choices = c("snpEff Annotated VCF (.vcf)", "Genotype Matrix (Processed)"),
            selected = character(0)
          ),
          uiOutput(ns("choice_output")),
        ),

        # Column Mapping Card
        bslib::card(
          bslib::card_header(tags$b("Column Mapping")),
          selectizeInput(ns("variant_id_col"),
            label = "Variant IDs Column",
            choices = NULL
          ),
          selectizeInput(ns("chrom_col"),
            label = "Chromosome Column",
            choices = NULL
          ),
          selectizeInput(ns("pos_col"),
            label = "Position Column",
            choices = NULL
          ),
          selectizeInput(ns("ref_al_col"),
            label = "Reference Allele Column",
            choices = NULL
          ),
          selectizeInput(ns("alt_al_col"),
            label = "Alternate Allele Column",
            choices = NULL
          ),
          numericInput(ns("geno_start"),
            label = "Genotype Data Start Column",
            value = 10
          )
        ),

        # Marker Selection Card
        bslib::card(
          bslib::card_header(tags$b("Marker Selection")),
          selectizeInput(ns("chr_ID"),
            label = "Chromosome ID",
            choices = NULL
          ),
          selectizeInput(ns("marker_ID"),
            label = "Marker ID",
            choices = NULL,
            multiple = TRUE,
            options = list(placeholder = "Select variants...")
          ),
          textInput(ns("reg_name"),
            label = "Region Name",
            placeholder = "e.g., drought resistance locus"
          )
        ),

        # Analysis Parameters Card
        bslib::card(
          bslib::card_header(tags$b("Analysis Parameters")),
          numericInput(ns("maf"),
            label = "Minor Allele Frequency (MAF)",
            value = 0.05, min = 0, max = 1
          ),
          bslib::input_switch(ns("draw_plot"),
            label = "Generate Alignment Plot",
            value = TRUE
          )
        ),
        bslib::card_footer(
          div(
            style = "text-align: center; margin-top: 15px;",
            actionButton(ns("run_but"),
              label = "Design Marker",
              icon = icon("drafting-compass"),
              class = "btn-primary",
              width = "100%"
            )
          )
        )
      ),
      # Main panel
      # Results container
      bslib::accordion(
        style = "margin-bottom: 70px;",
        id = ns("results_accordion"),
        width = "100%",
        open = TRUE,
        bslib::accordion_panel(
          "Comprehensive Table of KASP Marker Design Data and DNA Sequence Alignment to the Reference Genome",
          # preview for other generated markers
          selectizeInput(
            inputId = ns("done_markers"),
            label = "Select Marker ID",
            choices = NULL,
            width = "30%"
          ),
          DT::DTOutput(ns("kasp_table"), height = "200px"), # data.table output
          bslib::card(
            bslib::card_footer(
              fluidRow(
                column(width = 3, selectInput(
                  inputId = ns("exten"),
                  label = "Download file as?",
                  choices = c(".csv", ".xlsx"),
                  selected = ".xlsx",
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
        )
      ),

      # Conditional plot panel
      uiOutput(ns("plot_container"))
    )
  )
}

#' KASP Marker Design Server Function
#'
#' @description Server logic for KASP marker design module
#'
#' @param id Internal parameter for {shiny}.
#' @param input,output,session Internal parameters for {shiny}.
#' @importFrom shiny NS moduleServer observeEvent renderUI req reactiveVal
#' @importFrom shiny div icon tags strong
#' @noRd
#'
mod_kasp_marker_design_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


    # Server logic goes here
    observeEvent(input$upload_choice, {
      if (input$upload_choice == "snpEff Annotated VCF (.vcf)") {
        output$choice_output <- renderUI({
          fileInput(ns("vcf_file"),
            label = "Upload snpEff Annotated VCF (.vcf)",
            accept = ".vcf"
          )
        })
      } else {
        output$choice_output <- renderUI({
          fileInput(ns("gt_df"),
            label = "Upload Genotype Matrix (Processed)",
            accept = ".xlsx"
          )
        })
      }
    })


    # Read VCF data
    vcf_data <- reactive({
      if (is.null(input$vcf_file)) {
        return(NULL)
      }
      result <- marker.chr_ID(input$vcf_file$datapath)
      return(result)
    })

    # Update other parameter column for user.
    read_vcf_as_df <- function(vcf_file) {
      # Read VCF header to get column names
      header_lines <- readLines(vcf_file)
      header <- header_lines[grep("^#CHROM", header_lines)]
      colnames <- strsplit(header, "\t")[[1]]
      colnames[1] <- "CHROM" # fix formatting

      # Read VCF data (skip header lines)
      vcf_df <- utils::read.table(vcf_file,
        comment.char = "#", header = FALSE, sep = "\t",
        col.names = colnames, stringsAsFactors = FALSE
      )

      return(vcf_df)
    }

    # Reactive expression to populate inputs
    vcf_colnames <- reactive({
      req(input$vcf_file)
      colnames(read_vcf_as_df(input$vcf_file$datapath))
    })

    # Read Excel file
    gt_data <- reactive({
      req(input$gt_df)
      tryCatch(
        {
          result <- readxl::read_excel(input$gt_df$datapath)

          return(as.data.frame(result))
        },
        error = function(e) {
          shinyWidgets::show_toast(
            title = paste("Error reading Excel file:", e$message),
            type = "error",
            timer = 2000,
            position = "bottom-end"
          )
          NULL
        }
      )
    })


    # Extract column names
    gt_colnames <- reactive({
      req(gt_data())
      colnames(gt_data())[1:6]
    })

    # Extract unique chromosomes
    unique_chrom <- reactive({
      req(gt_data(), gt_colnames())

      unique(gt_data()[[gt_colnames()[2]]])
    })

    # Extract unique marker IDs
    unique_marker_id <- reactive({
      req(gt_data(), gt_colnames())

      unique(gt_data()[[gt_colnames()[1]]])
    })

    # Update column selection dropdowns
    observe({
      if (!is.null(input$vcf_file)) {
        # For VCF files

        updateSelectizeInput(
          server = TRUE, session, "variant_id_col",
          choices = vcf_colnames(),
          selected = vcf_colnames()[grep("id",
            x = vcf_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "chrom_col",
          choices = vcf_colnames(),
          selected = vcf_colnames()[grep("chro",
            x = vcf_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "pos_col",
          choices = vcf_colnames(),
          selected = vcf_colnames()[grep("pos",
            x = vcf_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "ref_al_col",
          choices = vcf_colnames(),
          selected = vcf_colnames()[grep("ref",
            x = vcf_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "alt_al_col",
          choices = vcf_colnames(),
          selected = vcf_colnames()[grep("alt",
            x = vcf_colnames(),
            ignore.case = TRUE
          )[1]]
        )
      } else if (!is.null(gt_data())) {
        # For Excel files
        updateSelectizeInput(
          server = TRUE, session, "variant_id_col",
          choices = gt_colnames(),
          selected = gt_colnames()[grep("id",
            x = gt_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "chrom_col",
          choices = gt_colnames(),
          selected = gt_colnames()[grep("chro",
            x = gt_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "pos_col",
          choices = gt_colnames(),
          selected = gt_colnames()[grep("pos",
            x = gt_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "ref_al_col",
          choices = gt_colnames(),
          selected = gt_colnames()[grep("ref",
            x = gt_colnames(),
            ignore.case = TRUE
          )[1]]
        )

        updateSelectizeInput(
          server = TRUE, session, "alt_al_col",
          choices = gt_colnames(),
          selected = gt_colnames()[grep("alt",
            x = gt_colnames(),
            ignore.case = TRUE
          )[1]]
        )
      }
    })

    # Update marker/chromosome dropdowns vcf
    observe({
      req(vcf_data())
      updateSelectizeInput(
        server = TRUE, session, "marker_ID",
        choices = vcf_data()$vcf_matrix_markerID
      )

      updateSelectizeInput(
        server = TRUE, session, "chr_ID",
        choices = vcf_data()$vcf_matrix_chromID
      )
    })

    # Update marker/chromosome dropdowns excel
    observe({
      req(gt_data())
      # Excel file logic
      req(unique_chrom(), unique_marker_id())

      updateSelectizeInput(
        server = TRUE, session, "marker_ID",
        choices = unique_marker_id()
      )
      updateSelectizeInput(
        server = TRUE, session, "chr_ID",
        choices = unique_chrom()
      )
    })

    # Reactive values
    kasp_des.result <- reactiveVal(NULL) # for dataframes
    kasp_des.plot <- reactiveVal(NULL) # for plots


    # Check which input is available and make use of it
    observeEvent(input$run_but, {
      list_markers <- list() # initialize list to store  marker dataframes
      list_plots <- list() # to store plot data frames

      req(
        input$variant_id_col, input$chrom_col,
        input$pos_col, input$ref_al_col, input$alt_al_col, input$geno_start,
        input$marker_ID, input$chr_ID, input$genome_file$datapath, input$maf
      )

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Designing KASP Marker... Please wait."
      )

      # First tryCatch for VCF file processing
      if (!is.null(input$vcf_file)) {
        tryCatch({
          for (marker in input$marker_ID) {
            # Using the altered kasp marker design here
            result_data <- kasp_marker_design(
              vcf_file = input$vcf_file$datapath,
              gt_df = NULL,
              variant_id_col = input$variant_id_col,
              chrom_col = input$chrom_col,
              pos_col = input$pos_col,
              ref_al_col = input$ref_al_col,
              alt_al_col = input$alt_al_col,
              geno_start = input$geno_start,
              marker_ID = marker,
              chr = input$chr_ID,
              genome_file = input$genome_file$datapath,
              plot_file = tempdir(),
              save_alignment = FALSE,
              region_name = input$reg_name,
              maf = input$maf
            )

            # store in list object created.
            list_markers[[marker]] <- result_data$marker_data
            list_plots[[marker]] <- result_data$plot
          }


          kasp_des.result(list_markers) # store dataframes
          kasp_des.plot(list_plots) # store  plots

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
        }, error = function(e) {
          kasp_des.result(NULL)
          shinyWidgets::show_alert(
            title = "Error!",
            text = paste("Failed to design KASP from VCF file:", e$message),
            type = "error",
            showCloseButton = TRUE,
            timer = 5000
          )
        }, finally = {
          shinybusy::remove_modal_spinner()
        })

        # Second tryCatch for gt_data processing
      } else if (!is.null(gt_data())) {
        tryCatch({
          for (marker in input$marker_ID) {
            result_data <- kasp_marker_design(
              vcf_file = NULL,
              gt_df = gt_data(),
              variant_id_col = input$variant_id_col,
              chrom_col = input$chrom_col,
              pos_col = input$pos_col,
              ref_al_col = input$ref_al_col,
              alt_al_col = input$alt_al_col,
              geno_start = input$geno_start,
              marker_ID = marker,
              chr = input$chr_ID,
              genome_file = input$genome_file$datapath,
              plot_file = tempdir(),
              region_name = input$reg_name,
              maf = input$maf,
              save_alignment = FALSE
            )

            list_markers[[marker]] <- result_data$marker_data
            list_plots[[marker]] <- result_data$plot
          }

          kasp_des.result(list_markers) # store dataframes
          kasp_des.plot(list_plots) # store plots

          # Success alert for gt_data processing
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
        }, error = function(e) {
          kasp_des.result(NULL)
          shinyWidgets::show_alert(
            title = "Error!",
            text = paste(
              "Failed to design KASP from genotype data:",
              e$message
            ),
            type = "error",
            showCloseButton = TRUE,
            timer = 5000
          )
        }, finally = {
          shinybusy::remove_modal_spinner()
        })
      }
    })

    # Update select input.
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

    # Render KASP table
    output$kasp_table <- DT::renderDT({
      req(kasp_des.result(), input$done_markers)
      DT::datatable(kasp_des.result()[[input$done_markers]],
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })


    # Download Table
    output$download_table <- downloadHandler(
      filename = function() {
        # Clean file name
        clean_name <- gsub("[^[:alnum:]_-]", "_", input$file_name)
        paste0(clean_name, input$exten)
      },
      content = function(file) {
        if (input$exten == ".csv") {
          write.csv(data.table::rbindlist(kasp_des.result()), file, row.names = FALSE)
        } else if (input$exten == ".xlsx") {
          openxlsx::write.xlsx(data.table::rbindlist(kasp_des.result()), file)
        }
      }
    )



    # Plot container UI - show if user selects true
    observeEvent(input$draw_plot, {
      if (input$draw_plot == TRUE) {
        output$plot_container <- renderUI({
          bslib::accordion(
            width = "100%",
            open = TRUE,
            bslib::accordion_panel(
              "KASP Sequence Alignment: 100 bp Upstream and Downstream of Target Site",
              selectizeInput(
                inputId = ns("plot_choice"),
                label = "Select Marker ID",
                width = "30%",
                choices = NULL
              ), # drop down for plots
              plotOutput(ns("plot"), height = "400px"),
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
# mod_kasp_marker_design_ui("kasp_marker_design_1")

## To be copied in the server
# mod_kasp_marker_design_server("kasp_marker_design_1")
