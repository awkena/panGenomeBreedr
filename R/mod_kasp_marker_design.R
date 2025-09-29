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
#' @importFrom bslib layout_sidebar sidebar card card_header card_footer input_switch accordion accordion_panel layout_columns
#'
mod_kasp_marker_design_ui <- function(id) {
  ns <- NS(id)
  tagList(
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 380,
        class = "bg-light",

        # Progress indicator
        div(
          class = "mb-3 p-2 rounded",
          style = "background-color: white; border-left: 4px solid #3498DB;",
          div(
            class = "d-flex align-items-center",
            icon("list-check", class = "text-info me-2"),
            strong("Design Workflow")
          )
        ),

        # Step 1: File Uploads
        bslib::accordion(
          id = "config_accordion",
          open = c("files", "mapping"),

          bslib::accordion_panel(
            title = div(
              icon("upload", class = "me-2"),
              "Step 1: Upload Files"
            ),
            value = "files",

            fileInput(ns("genome_file"),
                      label =   div(
                        icon("file-code", class = "me-2 text-success"),
                        "Genome Reference File"
                      ),
                      accept = c(".fa", ".fasta", ".gz")
            ),
            radioButtons(
              inputId = ns("upload_choice"),
              label = "Variant Data Type",
              choices = c("snpEff Annotated VCF (.vcf)", "Genotype Matrix (Processed)"),
              selected = character(0)
            ),
            uiOutput(ns("choice_output"))
          ),

          # Step 2: Column Mapping
          bslib::accordion_panel(
            title = div(
              icon("table-columns", class = "me-2"),
              "Step 2: Column Mapping"
            ),
            value = "mapping",

            div(
              class = "small text-muted mb-3",
              icon("info-circle", class = "me-1"),
              "Map your data columns to the required fields"
            ),

            selectizeInput(ns("variant_id_col"),
                           label = "Variant IDs Column",
                           choices = NULL
            ),

            bslib::layout_columns(
              col_widths = c(6, 6),
              selectizeInput(ns("chrom_col"),
                             label = "Chromosome:",
                             choices = NULL
              ),
              selectizeInput(ns("pos_col"),
                             label = "Position:",
                             choices = NULL
              )
            ),

            bslib::layout_columns(
              col_widths = c(6, 6),
              selectizeInput(ns("ref_al_col"),
                             label = "Reference Allele:",
                             choices = NULL
              ),
              selectizeInput(ns("alt_al_col"),
                             label = "Alternate Allele:",
                             choices = NULL
              )
            ),

            numericInput(ns("geno_start"),
                         label = "Genotype Data Start:",
                         value = 10
            )
          ),

          # Step 3: Marker Selection
          bslib::accordion_panel(
            title = div(
              icon("bullseye", class = "me-2"),
              "Step 3: Select Markers"
            ),
            value = "markers",

            selectizeInput(ns("chr_ID"),
                           label = "Target Chromosome:",
                           choices = NULL,
                           options = list(placeholder = "Choose chromosome...")
            ),

            selectizeInput(ns("marker_ID"),
                           label = "Marker Variants:",
                           choices = NULL,
                           multiple = TRUE,
                           options = list(placeholder = "Select one or more variants...")
            ),

            textInput(ns("reg_name"),
                      label = "Region Name:",
                      placeholder = "LGS1"
            ),

            div(
              class = "alert alert-info small mt-2",
              icon("lightbulb", class = "me-1"),
              "Give your region a descriptive name for easy identification"
            )
          ),

          # Step 4: Parameters
          bslib::accordion_panel(
            title = div(
              icon("sliders", class = "me-2"),
              "Step 4: Analysis Parameters"
            ),
            value = "params",

            numericInput(ns("maf"),
                         label = "Minor Allele Frequency (MAF):",
                         value = 0.05, min = 0, max = 1
            ),
            bslib::input_switch(ns("draw_plot"),
                                label = "Generate Alignment Plot",
                                value = TRUE
            )

          )
        ),

        bslib::card_footer(
          div(
            class = "mt-4 d-grid gap-2",
            actionButton(ns("run_but"),
                         label = "Design KASP Marker",
                         icon = icon("play", class = "me-2"),
                         #class = "btn-info",
                         class = "btn-success btn-lg",
                         style = "font-weight: 600; box-shadow: 0 4px 6px rgba(0,0,0,0.1);"
            )
          )
        )

      ),
      #-----------------------------------
      # Main panel Section
      #-----------------------------------
     div(class = "100",
         navset_card_pill(
           id = ns("results_tabs"),
           full_screen = TRUE,
           # Results table tab
           bslib::nav_panel(
             title = "Marker Data",
             icon = icon("table"),
             value = "table_tab",

             bslib::card(
               #   class = "h-100",
               bslib::card_header(
                 class = "bg-light",
                 div(
                   class = "d-flex justify-content-between align-items-center",
                   div(
                     icon("table", class = "me-2"),
                     strong("KASP Marker Design Results & Sequence Alignment")
                   )
                 )
               ),
               bslib::card_body(
                 class = "p-0",
                 DT::DTOutput(ns("kasp_table"),height = '500px'), # data.table output
               ),
               bslib::card_footer(
                 class = "bg-light",
                 bslib::layout_columns(
                   col_widths = c(3, 4,2),
                   selectInput(
                     inputId = ns("exten"),
                     label = "Download file as?",
                     choices = c(".csv", ".xlsx"),
                     selected = ".xlsx",
                     multiple = FALSE
                   ),
                   textInput(
                     inputId = ns("file_name"),
                     label = "Enter File Prefix",
                     value = "Kasp M_D for Intertek"
                   ),
                   div(
                     style = "display: flex; align-items: center; height: 100%; margin-top: 12px;",
                     downloadButton(
                       ns("download_table"),
                       label = "Export",
                       class = "btn-success w-100",
                       icon = icon("download")
                     )
                   )
                 )
               )
             )
           ),

           # Visualization tab
           bslib::nav_panel(
             title = "Alignment Plot",
             icon = icon("chart-bar"),
             value = "plot_tab",

             bslib::card(
               class = "h-100",
               card_header(
                 class = "bg-light",
                 div(
                   icon("chart-line", class = "me-2"),
                   strong("Sequence Alignment Visualization")
                 )
               ),
               bslib::card_body(
                 uiOutput(ns("plot_container"))
               ),
               bslib::card_footer(
                 class = "bg-light",
                 div(
                   class = "text-muted small",
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
mod_kasp_marker_design_server <- function(id){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

#---------
    #Dynamic rendering of file upload
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

#--------

    # Read VCF file and extract unique markers and chromosome ID
    vcf_data <- reactive({
       req(input$vcf_file)
       marker.chr_ID(input$vcf_file$datapath)

    })


    #Get column names of VCF file
    vcf_colnames <- reactive({
      req(input$vcf_file)
      colnames(read_vcf_as_df(input$vcf_file$datapath))
    })


#------------------------------
    #If user prefers to upload excel sheet containing genotype matrix
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
          return(NULL)
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

#----------------------

    # Reactively update select inputs depending on which file was available
    observeEvent(vcf_data() ,{
      req(vcf_data() ,  vcf_colnames())

        # For VCF files
        updateSelectizeInput(
          server = TRUE, session, "variant_id_col",
          choices = vcf_colnames(),
          selected = safe_grep_match(pattern = 'id',
                                     choices = vcf_colnames())[1]
        )

        updateSelectizeInput(
          server = TRUE, session, "chrom_col",
          choices = vcf_colnames(),
          selected = safe_grep_match(pattern = 'chro',
                                     choices = vcf_colnames())
        )

        updateSelectizeInput(
          server = TRUE, session, "pos_col",
          choices = vcf_colnames(),
          selected = safe_grep_match(pattern = 'pos',
                                     choices = vcf_colnames())
        )

        updateSelectizeInput(
          server = TRUE, session, "ref_al_col",
          choices = vcf_colnames(),
          selected = safe_grep_match(pattern = 'ref',
                                     choices = vcf_colnames())
        )

        updateSelectizeInput(
          server = TRUE, session, "alt_al_col",
          choices = vcf_colnames(),
          selected = safe_grep_match(pattern = 'alt',
                                    choices = vcf_colnames())
        )

        updateSelectizeInput(
          server = TRUE,
          session = session,
          inputId = "marker_ID",
          choices = vcf_data()$vcf_matrix_markerID
        )

        updateSelectizeInput(
          server = TRUE,
          session = session,
          inputId =  "chr_ID",
          choices = vcf_data()$vcf_matrix_chromID
        )

    })

    # Observer for matrix upload
    observeEvent(gt_data(),{
        req(unique_marker_id(),unique_chrom(), gt_colnames())
          # For Excel files
          updateSelectizeInput(
            server = TRUE, session, "variant_id_col",
            choices = gt_colnames(),
            selected = safe_grep_match(pattern = 'id',
                                       choices = gt_colnames())
          )

          updateSelectizeInput(
            server = TRUE, session, "chrom_col",
            choices = gt_colnames(),
            selected = safe_grep_match(pattern = "chro",
                                       choices = gt_colnames())
          )

          updateSelectizeInput(
            server = TRUE, session, "pos_col",
            choices = gt_colnames(),
            selected = safe_grep_match(pattern = 'pos',
                                       choices = gt_colnames())
          )

          updateSelectizeInput(
            server = TRUE, session, "ref_al_col",
            choices = gt_colnames(),
            selected = safe_grep_match(pattern = 'ref',
                                       choices = gt_colnames())
          )

          updateSelectizeInput(
            server = TRUE,
            session = session,
            inputId = "alt_al_col",
            choices = gt_colnames(),
            selected = safe_grep_match(pattern = 'alt',
                                       choices = gt_colnames())
          )

          updateSelectizeInput(
            server = TRUE,
            session =  session,
            inputId =  "marker_ID",
            choices = unique_marker_id()
          )
          updateSelectizeInput(
            server = TRUE,
            session =  session,
            inputId =  "chr_ID",
            choices = unique_chrom()
          )
      })

    # Empty reactive objects
    kasp_des.result <- reactiveVal(NULL) # Store dataframes
    kasp_des.plot <- reactiveVal(NULL) # store plots


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
        color =  "#27AE60",
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


          kasp_des.result(data.table::rbindlist(list_markers)) # store dataframes
          kasp_des.plot(list_plots) # store  plots

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

          kasp_des.result(data.table::rbindlist(list_markers))  # store dataframes
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


    # Update select input for user results preview
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
      req(kasp_des.result())
      DT::datatable(kasp_des.result(),
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })



    # Download table as batch
    output$download_table <- downloadHandler(
      filename = function() {
        # Clean file name
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
    observeEvent(input$draw_plot, {
      if (input$draw_plot == TRUE) {
        output$plot_container <- renderUI({
              tagList(
                selectizeInput(
                  inputId = ns("plot_choice"),
                  label = "Select Marker ID",
                  width = "30%",
                  choices = NULL
                ), # drop down for plots
                plotOutput(ns("plot"), height = "400px"),
                downloadButton(ns("download_plot"),
                               label = "Export All Plots (PDF)",
                               class = "btn-success w-30 mt-4",
                               width = "30%",
                               icon = icon("download")
                )
              )


        })
      } else if(input$draw_plot == FALSE) {
        output$plot_container <- renderUI({
          div(
            class = "text-center p-5 text-muted",
            icon("eye-slash", class = "fa-2x mb-3"),
            p("Plot generation disabled in parameters")
          )
        })
      }
    })

    # Render Plot
    output$plot <- renderPlot({

      req(kasp_des.plot(), input$plot_choice)
      print(kasp_des.plot()[[input$plot_choice]])
    })

    # Batch download plots
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
