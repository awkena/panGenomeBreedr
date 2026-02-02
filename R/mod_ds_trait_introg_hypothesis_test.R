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
  tagList(
    navset_card_underline(
      nav_panel(
        title = "Trait Introgression Hypothesis Testing",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(id = ns('sidebar'),
            width = 400,
            position = "left",
            class = "bg-light",

            # Accordion with organized sections
            bslib::accordion(
              id = "introgression_accordion",
              open = c("files", "mapfile", "settings", "qc"), # Panels open by default

              ## Data Acquisition
              bslib::accordion_panel(
                title = div(
                  icon("upload", class = "me-2"),
                  tags$span(
                    "Data Acquisition",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "files",

                fileInput(
                  inputId = ns("data_id"),
                  label = "Upload Genotype Data " ,
                  placeholder = ".csv,  .xlsx",
                  multiple = FALSE,
                  accept = c(".csv",".xlsx")
                ),

                layout_columns(col_widths = c(6,6),

                               selectInput(
                                 inputId = ns("data_type"),
                                 label = "Data Type:",
                                 choices = c("Agriplex", "Kasp / DArTag"),
                                 selected = "Agriplex",
                                 width = "100%"
                               ),
                               numericInput(inputId = ns("marker_start"),
                                            label = "Marker Start:",
                                            value = 7,
                                            min = 1,
                                            width = "100%" )
                               ),
                textInput(
                  inputId = ns("allele_sep"),
                  label = "Enter Allele Separator for Data Format",
                  value = " / ",
                  width = "100%"
                )
              ),

              # Genomic Mapping Setup
              bslib::accordion_panel(
                title = div(
                  icon("map", class = "me-2"),
                  tags$span(
                    "Genomic Mapping Setup",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "mapfile",

                radioButtons(
                  inputId = ns("choice"),
                  label = "Do you have a map file?",
                  choices = c("Yes" = "yes", "Generate one" = "no"),
                  selected = "yes"
                ),

                # Dynamic rendering based on choice
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
                      label = "Select Column for SNP ID",
                      choices = NULL,
                      width = "100%"
                    )
                  )
                ),

                conditionalPanel(
                  condition = paste0('input["', ns("choice"), '"] == "no"'),
                  tagList(
                    textInput(
                      inputId = ns("prefix_marker"),
                      label = "Enter Prefix for Marker ID",
                      value = "S",
                      width = "100%"
                    ),
                    textInput(
                      inputId = ns("sep_marker"),
                      label = "Enter Separator for Marker ID",
                      value = "_",
                      width = "100%"
                    )
                  )
                )
                # div(
                #   class = "alert alert-info small mt-2",
                #   icon("lightbulb", class = "me-1"),
                #   "Map file links markers to their genomic positions"
                # )
              ),

              ## Selection Parameters
              bslib::accordion_panel(
                title = div(
                  icon("sliders", class = "me-2"),
                  tags$span(
                    "Selection Parameters",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "settings",

                selectizeInput(
                  inputId = ns("sample_id"),
                  label = "Select Sample ID Column",
                  choices = NULL,
                  width = "100%"
                ),

                  selectizeInput(
                    inputId = ns("cluster_by"),
                    label = "Group By:",
                    choices = NULL,
                    width = "100%"
                  ),

                  uiOutput(outputId = ns("cluster_focus_ui")),

                bslib::layout_columns(
                  col_widths = c(6, 6),
                  div(
                    numericInput(
                      inputId = ns("rp"),
                      label = "Recurrent Parent Index:",
                      value = 1,
                      min = 1,
                      width = "100%"
                    ),
                    uiOutput(outputId = ns('recurrent_help'))
                  ),
                  div(
                    numericInput(
                      inputId = ns("dp"),
                      label = "Donor Parent Index:",
                      value = 3,
                      min = 1,
                      width = "100%"
                    ),
                    uiOutput(outputId = ns('donor_help'))
                  )

                )
              ),

              ## Genotypic Data Cleaning
              bslib::accordion_panel(
                title = div(
                  icon("filter", class = "me-2"),
                  tags$span(
                    "Genotypic Data Cleaning",
                    style = "font-weight: bold; font-size: 1.1rem;"
                  )
                ),
                value = "qc",

                div(
                  class = "small text-muted mb-3",
                  icon("info-circle", class = "me-1"),
                  "Enable filters to improve data quality"
                ),

                bslib::input_switch(
                  id = ns("apply_par_poly"),
                  label = "Remove Monomorphic Parents",
                  value = TRUE
                ),

                bslib::input_switch(
                  id = ns("apply_par_miss"),
                  label = "Remove Missing Parent Data",
                  value = TRUE
                ),

                bslib::input_switch(
                  id = ns("apply_geno_good"),
                  label = "Apply Genotype Error Check",
                  value = TRUE
                ),

                bslib::input_switch(
                  id = ns("apply_par_homo"),
                  label = "Filter Heterozygous Parents",
                  value = TRUE
                )
              )
            ),

            # Action button to run analysis
            div(
              class = "mt-4 d-grid gap-2",
              actionButton(
                inputId = ns("config"),
                label = "Process Data",
                icon = icon("gears", class = "me-2"),
                class = "btn-primary btn-lg",
                style = "font-weight: 600;"
              )
            )
          ),

          # Main content area
          bslib::navset_hidden(id = ns('config_pages'),
            bslib::nav_panel_hidden(
              value = "main_config",
              bslib::input_switch(id = ns("configure"), label = "Plot Configuration", value = T),
              conditionalPanel(
                condition = paste0('input["', ns("configure"), '"] == true'),
                card(
                  fluidRow(
                    # Card 1: Heatmap Configuration
                    column(
                      width = 4,
                      bslib::card(
                        class = "shadow p",
                        max_height = "850px",
                        height = "100%",
                        bslib::card_header(tags$b("Heatmap Computation Parameters"),
                                           class = "bg-warning text-center",
                                           style = "font-size:18px;"
                        ),
                        bslib::card_body(
                          selectInput(
                            inputId = ns("snp_ids"),
                            label = "Select Marker Column",
                            choices = NULL,
                            width = "100%"
                          ),
                          selectInput(
                            inputId = ns("chr"),
                            label = "Select Chromosome Column",
                            choices = NULL,
                            width = "100%"
                          ),
                          selectInput(
                            inputId = ns("chr_pos"),
                            label = "Select Position Column",
                            choices = NULL,
                            width = "100%"
                          ),
                          selectInput(
                            inputId = ns("parents"),
                            label = "Choose Parents",
                            choices = NULL,
                            width = "100%",
                            multiple = TRUE
                          ),
                          numericInput(
                            inputId = ns("group_sz"),
                            label = "Progeny Batch Size",
                            value = 0,
                            min = 0,
                            step = 1,
                            width = "100%"
                          ),
                          radioButtons(
                            inputId = ns("options"),
                            label = "Display Annotations?",
                            choices = c(
                              "Hide Annotations" = "no",
                              "Show Annotations" = "yes" # cross_qc_annotate()/ # cross_qc_heatmap2()
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
                              class = "btn-outline-primary",
                              style = "font-weight: bold;"
                            )
                          ),
                          radioButtons(
                            inputId = ns("focus_qtl"),
                            label = "Annotate Specific QTL?",
                            choices = c("Yes" = "yes", "No" = "no"),
                            selected = "no",
                            inline = TRUE
                          ),
                          conditionalPanel(
                            # reveal based on choice -- Annotate specific QTL
                            condition = "input.focus_qtl == 'yes'",
                            ns = ns,
                            div(
                              selectInput(
                                inputId = ns("qtls"),
                                label = "Select a QTL",
                                choices = NULL,
                                multiple = T
                              ),
                              # Using span to keep the icon and text on the same line
                              helpText(
                                span(icon("info-circle", vertical_align = "middle")),
                                "QTL coordinates must first be defined!"
                              )
                            )
                          )
                        )
                      )
                    ),

                    # Card 2: Heatmap Visualization Controls
                    column(
                      width = 8,
                      bslib::card(
                        class = "shadow p",
                        max_height = "850px",
                        height = "100%",
                        bslib::card_header(tags$b("Heatmap Visualization Controls"),
                                           class = "bg-info text-center",
                                           style = "font-size:18px;"
                        ),
                        bslib::card_body(
                          fluidRow(
                            # Left Column - 4 widgets
                            column(
                              width = 6,
                              textInput(
                                inputId = ns("legend_title"),
                                label = "Enter Plot Legend Title",
                                value = "Heatmap Key",
                                width = "100%"
                              ),
                              selectInput(
                                inputId = ns("col_mapping"),
                                label = "Choose Heatmap Colors",
                                choices = grDevices::colors(),
                                multiple = TRUE,
                                width = "100%"
                              ),
                              selectInput(
                                inputId = ns("col_labels"),
                                label = "Genotype Labels",
                                choices = NULL,
                                multiple = TRUE,
                                width = "100%"
                              ),
                              selectInput(
                                inputId = ns("panel_fill"),
                                label = "Panel Background Fill Color",
                                choices = grDevices::colors(),
                                selected = "grey80",
                                width = "100%"
                              ),
                              numericInput(
                                inputId = ns("text_scale_fct"),
                                label = "Text Scaling Size",
                                value = 0.3,
                                min = 0.1,
                                step = 0.01,
                                max = 1
                              )
                            ),
                            # Right Column - 4 widgets
                            column(
                              width = 6,
                              selectInput(
                                inputId = ns("panel_col"),
                                label = "Panel Border Color",
                                choices = grDevices::colors(),
                                selected = "white",
                                width = "100%"
                              ),
                              numericInput(
                                inputId = ns("alpha"),
                                label = "Point Transparency",
                                value = 1,
                                min = 0,
                                max = 1,
                                step = 0.1,
                                width = "100%"
                              ),
                              numericInput(
                                inputId = ns("text_size"),
                                label = "Text Size",
                                value = 12,
                                min = 1,
                                step = 0.1,
                                width = "100%"
                              ),
                              numericInput(
                                inputId = ns("label_offset"),
                                label = "Trait Label Position",
                                value = 0.4,
                                min = -10,
                                max = 10,
                                step = 0.1,
                                width = "100%"
                              )
                            )
                          )
                          # helpText(tags$b("Key for ordering heatmap colors")),
                          # verbatimTextOutput(outputId = ns("color_format"))
                        )
                      )
                    )
                  ),
                  actionButton(
                    inputId = ns("generate_heatmap"),
                    label = "Generate Heatmap",
                    icon = icon("play", class = "me-2"),
                    class = "btn-success btn-lg",
                    style = "font-weight: 600;"
                  )
                )
              )

            ),

            bslib::nav_panel_hidden(
              value = 'coord_entry',
              bslib::card(
                full_screen = TRUE,
                bslib::card_header(
                  div(class = "d-flex justify-content-between align-items-center",
                      tags$b("Genomic Coordinate Definition"),
                      actionButton(ns("addTrait"), "Add New Trait", icon = icon("plus"), class = "btn-success btn-sm"))
                ),
                bslib::card_body(
                  p(class = "text-muted", "Define specific points or genomic ranges for heatmap annotation."),
                  div(id = ns("inputs_container"),
                      makeTraitInput(1, ns) # Initial row
                  )
                ),
                bslib::card_footer(
                  div(class = "d-flex justify-content-end gap-2",
                      actionButton(ns("cancel_coords"), "Cancel", class = "btn-light"),
                      actionButton(ns("submit_coords"), "Save & Apply", icon = icon("check"), class = "btn-primary"))
                )
              )
            )
          ),

          # Results Section
          fluidRow(
            column(
              width = 12,
              bslib::accordion(height = "800px",
                bslib::accordion_panel(
                  title = "Heatmap Results & Analysis",
                  icon = icon("chart-line"),
                  div(
                    selectInput(
                      inputId = ns("batch_no"),
                      label = "Select Batch for Display",
                      choices = NULL,
                      width = "30%"
                    ),
                    plotOutput(
                      outputId = ns("ant_heatmap"),
                      width = "100%",
                      height = "650px"
                    ),
                    bslib::card(card_footer(
                      fluidRow(
                        column(
                          3,
                          textInput(
                            inputId = ns("file_name"),
                            label = "Enter Filename",
                            value = "Heatmap_result"
                          )
                        ),
                        column(
                          3,
                          numericInput(
                            inputId = ns("width"),
                            label = "Set Plot Width",
                            value = 10, min = 1
                          )
                        ), column(
                          3,
                          numericInput(
                            inputId = ns("height"),
                            label = "Set Plot Height",
                            value = 7, min = 1
                          )
                        )
                      ),
                      downloadButton(
                        outputId = ns("download_plot1"),
                        label = "Download Plot", class = "btn-success"
                      )
                    ))
                  )
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

    # # Dynamic ui based on choice of individual.
    observe({
      if (input$choice == "no") {
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

              # Using a list with custom spacing
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

              # Example to user, highlighted in a Card-style div
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




    # Read and validate uploaded CSV file
    validated_data <- eventReactive(input$data_id, {
      req(input$data_id)
       read_mapfile(filepath = input$data_id$datapath)
    })



    # Get column names from the validated data
    data_colnames <- reactive({
      req(validated_data())
      colnames(validated_data())
    })


    # Populate Batch and Genotype column selectors dynamically
    observe({
      req(data_colnames())

      # Update batch column selector
      updateSelectizeInput(
        session,
        inputId = "cluster_by",
        choices = c('None',data_colnames()),server = T
      )

      # Update genotype column selector
      updateSelectizeInput(
        session,
        inputId = "sample_id",
        choices = data_colnames(),
        selected = safe_grep_match(
          pattern = "genotype",
          choices = data_colnames()
        ),server = T
      )
    })


    # cluster logic
    output$cluster_focus_ui <- renderUI({
      # Only show the second dropdown if a factor is selected
      req(input$cluster_by , validated_data())
      if (input$cluster_by == "None") return(NULL)

      # Get the unique values
      choices_to_focus <- unique(validated_data()[[input$cluster_by]])

      selectInput(inputId = ns("cluster_focus"),
                  label = paste("Select", input$cluster_by, "to focus on:"),
                  choices = choices_to_focus)
    })




     # Update the validated data based on the cluster selected.
    sub_validated_data <- reactive({
      req(validated_data(),input$cluster_by)

      # Check to prevent undefined columns from being selected.
      if (!input$cluster_by %in% colnames(validated_data())){
        return(validated_data())
      } else if(input$cluster_by != 'None'){
        req(input$cluster_focus)
        validated_data()[validated_data()[input$cluster_by] == input$cluster_focus ,]
      }
    })


    ## Help text for donor and reccurent parents
    # Display the Donor Name
      output$donor_help <- renderUI({
        if(input$cluster_by != "None"){
        req(input$dp ,input$sample_id , sub_validated_data())
        name <- sub_validated_data()[[input$sample_id]][input$dp]
        helpText(strong("Selected: "), span(name, style = "color: #e67e22;"))
        }
        else{
          # if user's data cannot be partitioned into groups
          req(validated_data(), input$sample_id,input$dp )
          name <- validated_data()[input$dp,input$sample_id]
          helpText(strong("Selected: "), span(name, style = "color: #0275d8; font-size: 1.1em;"))
        }
      })

      #  Display the Recurrent Name
      output$recurrent_help <- renderUI({
        if(input$cluster_by != "None"){
          req(input$rp ,input$sample_id , sub_validated_data())
          name <- sub_validated_data()[[input$sample_id]][input$rp]
          helpText(strong("Selected: "), span(name, style = "color: #0275d8; font-size: 1.1em;"))
        } else{
          # if user's data cannot be partitioned into groups
          req(validated_data(), input$sample_id,input$rp )
          name <- validated_data()[input$rp,input$sample_id]
          helpText(strong("Selected: "), span(name, style = "color: #0275d8; font-size: 1.1em;"))
        }

      })



    # Read map file if user has.
    map_file <- reactive({
      req(input$mapfile)
      read_mapfile(filepath = input$mapfile$datapath)
    })

    # Populate field for snp id column based on mapfile
    observe({
      req(map_file())
      updateSelectInput(session,
                        inputId = "snp_id",
                        choices = colnames(map_file()),
                        selected = safe_grep_match(pattern = 'id',
                                                   choices = colnames(map_file()))
      )
    })


    # process data.
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
                         rp = input$rp,
                         dp = input$dp,
                         Prefix = if (input$choice == "no") input$prefix_marker else NULL,
                         feedback = input$choice,
                         mapfile_path = if (!is.null(map_file())) map_file() else NULL)

      }, error = function(e) {
        shinyWidgets::show_alert(
          title = 'Error',
          text = paste("An error occurred while processing the data:", e$message),
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

   # Update user on next step
   observeEvent(Result(), {
     req(Result()) # Only proceed if Result is not NULL

     #bslib::sidebar_toggle(id = 'sidebar' , open = "closed")
     shinyWidgets::show_alert(
       title = "Data Processed!",
       text = "Your genotype data is ready. Please proceed to 'Plot Configuration' to generate your heatmap.",
       type = "success",
       btn_labels = "Got it!",
       closeOnClickOutside = FALSE,
       showCloseButton = TRUE
     )
   })



   # Obtain the order of labels
    color_code_checker_res <- reactive({
      req(Result() , input$dp,input$rp )
      color_code_checker(Result()$proc_kasp_f ,
                         par_1 = rownames(Result()$proc_kasp_f)[input$rp],
                         par_2 = rownames(Result()$proc_kasp_f)[input$dp])
    })



    # Observer for updating map-related inputs
    observe({
      req(Result())

      # Update SNP ID selection
      updateSelectInput(session,
                        inputId = "snp_ids",
                        choices = colnames(Result()$mapfile),
                        selected = safe_grep_match(pattern = 'id',
                                                   choices = colnames(Result()$mapfile))

      )

      # Update Chromosome selection
      updateSelectInput(session,
                        inputId = "chr",
                        choices = colnames(Result()$mapfile),
                        selected = safe_grep_match(pattern = 'chr',
                                                   choices = colnames(Result()$mapfile))
      )

      # Update Position selection
      updateSelectInput(session,
                        inputId = "chr_pos",
                        choices = colnames(Result()$mapfile),
                        selected = safe_grep_match(pattern = 'pos',
                                                   choices = colnames(Result()$mapfile))
      )

      req(color_code_checker_res())
      # Update Position selection
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


    ### Heat MAP Plot Section  ###
    values <- reactiveValues(
      count = 1,
      removed_traits = c(),
      trait_data = NULL
    )

    # Navigation Logic
    observeEvent(input$newBtn, {
      bslib::nav_select("config_pages", "coord_entry")
    })

    # Adding QTLs
    observeEvent(input$addTrait, {
      values$count <- values$count + 1
      insertUI(
        selector = paste0("#", ns("inputs_container")),
        where = "beforeEnd",
        ui = makeTraitInput(values$count, ns)
      )
    })

    # Removing added QTL
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

    # Submit
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
        shinyWidgets::show_alert("Error", "Please define at least one valid trait with a name and type.", "error")
        return()
      }

      values$trait_data <- trait_list
      bslib::nav_select("config_pages", "main_config")

      shinyWidgets::show_toast("Success", "Trait coordinates saved.", type = "success")
    })

    # Cancel button
    observeEvent(input$cancel_coords, {
      bslib::nav_select("config_pages", "main_config")

      values$trait_data <- NULL # nullify all existing trait data

      # Toast to confirm nothing was saved
      shinyWidgets::show_toast(
        title = "Cancelled",
        text = "No changes were made to trait coordinates.",
        type = "info"
      )
    })


    ## User defined qtl
    observe({
      req(values$trait_data)
      if(input$focus_qtl == 'yes'){
        updateSelectInput(session, inputId = 'qtls',choices = names(values$trait_data))
      } else {
        updateSelectInput(session, inputId = 'qtls',choices = NULL,selected = NULL)
      }
    })


    ## Heat map generation.
    heatmap_plot <- eventReactive( input$generate_heatmap,{
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

    # Update progeny batch size drop down..
    observe({
      req(heatmap_plot())
      updateSelectInput(session,inputId = 'batch_no',choices = names(heatmap_plot()))
    })



    # Update map dynamically.
    heatmap_1 <- reactive({
      req(heatmap_plot(), input$batch_no)
      if (input$batch_no %in% names(heatmap_plot())) {
        heatmap_plot()[[input$batch_no]]
      } else {
        NULL
      }
    })

    # Print plot
    output$ant_heatmap <- renderPlot({
      tryCatch({
        plot_obj <- heatmap_1()
        req(plot_obj)

        if (inherits(plot_obj, "ggplot")) {
          print(plot_obj)
          # show toast
          shinyWidgets::show_toast(
            title = '',
            text = 'Heatmap Generated Successfully',
            type = "success"
          )

        } else {

          plot_obj
        }

      }, error = function(e) {
        return(NULL)
        shinyWidgets::show_alert(
          title = 'Error',
          text = paste("An error occurred while generating plot", e$message),
          type = "error"
        )
      })
    })



    # Download batch plot.
    output$download_plot1 <- downloadHandler(
      filename = function() {
        req(input$file_name)
        paste0(input$file_name, ".pdf")  # Single PDF containing ALL plots
      },
      content = function(file) {
        req(heatmap_plot())  # Ensure plots exist

        tryCatch(
          {
            # Start PDF (onefile=TRUE ensures multi-page)
            grDevices::pdf(file,
                           width = input$width,
                           height = input$height,
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


