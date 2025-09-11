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
          sidebar = bslib::sidebar(
            width = 400,
            position = "left",
            bslib::card(
              height = "100%",
              bslib::card_body(
              bslib::card(
                bslib::card_header(tags$b('Upload Input Files'),
                                   class = 'text-center',
                                   style = "font-size:18px;"),
                bslib::card_body(
                  fileInput(
                    inputId = ns("data_id"),
                    label = "Upload Kasp/Agriplex File",
                    multiple = FALSE,
                    accept = ".csv",
                    width = "100%"
                  ),
                  selectInput(
                    inputId = ns("data_type"),
                    label = "Indicate Data Format",
                    choices = c("agriplex", "Kasp"),
                    selected = "agriplex",
                    width = "100%"
                  ), # user defines data format

                  textInput(
                    inputId = ns("allele_sep"),
                    label = "Enter Allele Separator for Data Format",
                    value = " / ",
                    width = "100%"
                  )
                )
              ),
              bslib::card(
                bslib::card_header(tags$b('Mapfile Setup'),
                                   class = 'text-center',
                                   style = "font-size:18px;"),
                bslib::card_body(
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
                        accept = c(".csv",".xlsx",".xls"),
                        width = "100%"
                      ),
                      selectInput(
                        inputId = ns("snp_id"), #
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
                )
              ),

              bslib::card(
                bslib::card_header(tags$b('Genotype & Batch Settings'),
                                   class = 'text-center',
                                   style = "font-size:18px;"),
                bslib::card_body(
                  # Batch  & genotype settings.
                  selectInput(
                    inputId = ns("genotype_col"),
                    label = "Select Genotype Column",
                    choices = NULL,
                    width = "100%"
                  ),
                  selectInput(
                    inputId = ns("batch_col"),
                    label = "Select Batch Column",
                    choices = NULL,
                    width = "100%"
                  ), # select batch column will populates

                  selectInput(
                    inputId = ns("batch"),
                    label = "Select Focused Batch",
                    choices = NULL,
                    width = "100%"
                  ), # unique batches will come here.

                  uiOutput(ns("marker_sep")), # this must be dynamic based on users choice

                  selectInput(
                    inputId = ns("dp"),
                    label = "Select Donor Parent",
                    choices = NULL,
                    width = "100%"
                  ),
                  selectInput(
                    inputId = ns("rp"),
                    label = "Select Recurrent Parent",
                    choices = NULL,
                    width = "100%"
                  )
                )
                )
              ),
              bslib::card(
                bslib::card_header(tags$b('Quality Control Switches'),
                                   class = 'text-center',
                                   style = "font-size:18px;"),
                bslib::card_body(
                  tagList(
                    bslib::input_switch(
                      id = ns("apply_par_poly"),
                      label = "Remove Monomorphic Parents",
                      value = TRUE,
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
                )
              ),
              card_footer(
                actionButton(
                  inputId = ns("config"),
                  label = "Submit",
                  icon = icon("check"),
                  width = "100%",
                  #class = "btn-info"
                  style = "background-color: forestgreen; color: white; font-weight: bold; border: none;",
                  `onmouseover` = "this.style.backgroundColor='#145214'",
                  `onmouseout` = "this.style.backgroundColor='forestgreen'"
                )
              )
            )
          ),

          # Main content area
          bslib::input_switch(id = ns("configure"), label = "Configure Plot", value = FALSE),
          conditionalPanel(
            condition = paste0('input["', ns("configure"), '"] == true'),
            fluidRow(
              # Card 1: Heatmap Configuration
              column(
                width = 4,
                bslib::card(
                  class = "shadow p",
                  max_height = "850px",
                  height = "100%",
                  bslib::card_header(tags$b("Heatmap Configuration"),
                                     class = "bg-success text-center",
                                     style = "font-size:18px;"
                                     ),
                  bslib::card_body(
                    selectInput(
                      inputId = ns("snp_ids"),
                      label = "Select Column for SNP ID",
                      choices = NULL,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("chr"),
                      label = "Select Column for Chromosome",
                      choices = NULL,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("chr_pos"),
                      label = "Select Column for Chromosome Position",
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
                      label = "Batch Size of Progenies to Include in Heatmap",
                      value = 0,
                      min = 0,
                      step = 1,
                      width = "100%"
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
                          value = 0.18,
                          min = 0.1,
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
                          width = "100%"
                        ),
                        numericInput(
                          inputId = ns("label_offset"),
                          label = "Trait Label Position",
                          value = -1,
                          min = -10,
                          width = "100%"
                        )
                      )
                    ),
                    helpText(tags$b('Key for ordering heatmap colors')),
                    verbatimTextOutput(outputId = ns('color_format'))

                  )
                )
              )
            )
          ),

          # Results Section
          fluidRow(
            column(
              width = 12,
              bslib::accordion(
                bslib::accordion_panel(
                  title = "Heatmap Results & Analysis",
                  icon = icon("chart-line"),
                  bslib::navset_card_tab(
                    bslib::nav_panel(
                      title = "Annotated Heatmaps",
                      icon = icon("tags"),
                      radioButtons(
                        inputId = ns("choices"),
                        label = "Select Choice",
                        choices = c(
                          "Default" = "def",
                          "Window Size" = "widsz",
                          "Annotate" = "antt"
                        ),
                        selected = "def",
                        inline = TRUE,
                        width = "500px"
                      ),
                      # Conditional panels based on radio button selection
                      conditionalPanel(
                        condition = paste0("input['", ns("choices"), "'] == 'def'"),
                        nav_panel(
                          title = "Default",
                          selectInput(
                            inputId = ns("batch_no"),
                            label = "Select Batch for Display",
                            choices = NULL,
                            width = "30%"
                          ),
                          plotOutput(
                            outputId = ns("heatmap"),
                            width = "100%",
                            height = "600px"
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
                                  value = 14, min = 1
                                )
                              ), column(
                                3,
                                numericInput(
                                  inputId = ns("height"),
                                  label = "Set Plot Height",
                                  value = 6, min = 1
                                )
                              )
                            ),
                            downloadButton(
                              outputId = ns("download_plot1"),
                              label = "Download Plot", class = "btn-success"
                            )
                          ))
                        )
                      ),

                      conditionalPanel(
                        condition = paste0("input['", ns("choices"), "'] == 'widsz'"),
                        nav_panel(
                          title = "Window Size",
                          selectInput(
                            inputId = ns("chr_vw"),
                            label = "Select Chromosome to View",
                            choices = NULL,
                            width = "30%"
                          ),
                          plotOutput(
                            outputId = ns("wind_heatmap"),
                            width = "100%",
                            height = "600px"
                          )
                        )

                      ),

                      conditionalPanel(
                        condition = paste0("input['", ns("choices"), "'] == 'antt'"),
                        nav_panel(
                          title = "Annotation",
                          div(
                            style = "display: flex; align-items: center; gap: 20px; flex-wrap: wrap;",
                            actionButton(
                              inputId = ns("newBtn"),
                              label = "Set Trait Positions",
                              icon = icon("plus-circle"),
                              style = "background-color: forestgreen; color: white; font-weight: bold; border: none; flex-shrink: 0;",
                              `onmouseover` = "this.style.backgroundColor='#145214'",
                              `onmouseout` = "this.style.backgroundColor='forestgreen'"
                            ),
                            radioButtons(
                              inputId = ns("options"),
                              label = "",
                              choices = c(
                                "Box Annotation" = "1",
                                "Line Annotation" = "2"
                              ),
                              selected = "1",
                              inline = TRUE
                            )
                          ),
                          plotOutput(
                            outputId = ns("ant_heatmap"),
                            width = "100%",
                            height = "700px"
                          ),
                          bslib::card(card_footer(
                            fluidRow(
                              column(
                                3,
                                textInput(
                                  inputId = ns("file_name2"),
                                  label = "Enter Filename",
                                  value = "Heatmap_result"
                                )
                              ),
                              column(
                                3,
                                numericInput(
                                  inputId = ns("width2"),
                                  label = "Set Plot Width",
                                  value = 8, min = 1
                                )
                              ), column(
                                3,
                                numericInput(
                                  inputId = ns("height2"),
                                  label = "Set Plot Height",
                                  value = 6, min = 1
                                )
                              )
                            ),
                            downloadButton(
                              outputId = ns("download_plot2"),
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
      ),
      nav_panel(
        title = "Data QC", icon = icon("bell"),
        splitLayout(
          bslib::card(
            card_header(tags$b("SNP Loci with Potential Genotype Call Errors")),
            bslib::card_body(
              DT::DTOutput(outputId = ns("geno_error_tbl"))
            )
          ),
          bslib::card(
            card_header(tags$b("SNP Loci with Parent Missing")),
            bslib::card_body(
              DT::DTOutput(outputId = ns("par_missing_tbl"))
            )
          )
        ),
        splitLayout(
          bslib::card(
            card_header(tags$b("SNP Loci with Heterozygote Parent")),
            bslib::card_body(
              DT::DTOutput(outputId = ns("par_het_tbl"))
            )
          ),
          bslib::card(
            card_header(tags$b("SNP Loci with Unexpected Locus")),
            bslib::card_body(
              DT::DTOutput(outputId = ns("unexp_locu_tbl"))
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
            title = tags$b("Important Note!!"),
            tagList(
              p("Marker names must follow a structured format to be parsed into the map file."),
              tags$ul(
                tags$li("A common prefix before each marker (e.g. ", tags$b("S"), ")."),
                tags$li("Chromosome number immediately after the prefix (e.g. ", tags$b("1"), ")."),
                tags$li("A separator character (e.g. ", tags$b("_"), ")."),
                tags$li("Position number after the separator (e.g. ", tags$b("101"), ").")
              )
              ,
              p(
                "Example: ", tags$b("S1_101"),
                ".Where ", tags$b("S"), " is the prefix, ",
                tags$b("1"), " is the chromosome number, and ",
                tags$b("101"), " is the position."
              )

            ),
            easyClose = FALSE,
            footer = modalButton("Got it!")
          )
        )
      }
    })

    # Read and validate uploaded CSV file
    validated_data <- eventReactive(input$data_id, {
      req(input$data_id)

      data <- read_mapfile(filepath = input$data_id$datapath)

      # Validate required column names
      tryCatch({
        check_colnames_validate(data)
        return(data)  # Return validated data
      }, error = function(e) {
        shinyWidgets::show_alert(
          title = "Column Validation Error",
          text = e$message,
          type = "error"
        )
        NULL
      })
    })



    # Get column names from the validated data
    data_colnames <- reactive({
      req(validated_data())
      Get_dt_coln(validated_data())
    })


    # Populate Batch and Genotype column selectors dynamically
    observe({
      req(data_colnames())

      # Update batch column selector
      updateSelectInput(
        session,
        inputId = "batch_col",
        choices = data_colnames(),
        selected = safe_grep_match(
          pattern = "batch",
          choices = data_colnames()
        )
      )

      # Update genotype column selector
      updateSelectInput(
        session,
        inputId = "genotype_col",
        choices = data_colnames(),
        selected = safe_grep_match(
          pattern = "genotype",
          choices = data_colnames()
        )
      )
    })


    # Get unique batch values from the validated dataset
    uniq_batch <- reactive({
      req(validated_data(), input$batch_col)
      unique(validated_data()[[input$batch_col]])
    })

    # Populate batch selector input
    observe({
      req(uniq_batch())
      updateSelectInput(session,
                        inputId = "batch",
                        choices = uniq_batch(),
                        selected = uniq_batch()[1]
      )
    })

    # Get genotypes and populate for parents.
    Genotype_names <- reactive({
      req(input$batch, validated_data())
      Genotypes_user(data = validated_data(), Batch = input$batch )
    })

    # Update recurrent and donor parent selector inputs
    observe({
      req(Genotype_names())
      updateSelectInput(session,
                        inputId = "dp",
                        choices = Genotype_names())

      updateSelectInput(session,
                        inputId = "rp",
                        choices = Genotype_names(),
                        selected = Genotype_names()[3])

      updateSelectInput(session,
                        inputId = "parents",
                        choices = Genotype_names(),
                        selected = Genotype_names()[c(1, 3)]
      )
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


    # store the parents selection temporary
    values <- reactiveValues(
      locked_parents = NULL
    )


    # process data.
   Result <- eventReactive(input$config, {
      req(
        validated_data(), input$batch, input$batch_col,
        input$data_type, input$allele_sep,
        input$rp, input$dp, input$choice, data_colnames(),
        input$genotype_col, Genotype_names()
      )

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#0dc5c1",
        text = "Getting results... Please wait."
      )

      result <- tryCatch({
        proc_nd_map_func(
          data = validated_data(),
          Batch = input$batch,
          batch_col = input$batch_col,
          marker_sep = if (input$choice == "no") input$sep_marker else NULL,
          apply_par_poly = input$apply_par_poly,
          apply_par_miss = input$apply_par_miss,
          apply_geno_good = input$apply_geno_good,
          apply_par_homo = input$apply_par_homo,
          genotype = input$genotype_col,
          snp_id = if (input$choice == "yes") input$snp_id else NULL,
          calls_sep = check_sep(input$allele_sep),
          data_type = input$data_type,
          rp = input$rp,
          dp = input$dp,
          Prefix = if (input$choice == "no") input$prefix_marker else NULL,
          geno_vec = Genotype_names(),
          feedback = input$choice,
          na_code = NA,
          data_col = data_colnames(),
          mapfile_path = if (!is.null(map_file())) map_file() else NULL
        )

      }, error = function(e) {
        shinyWidgets::show_alert(
          title = 'Error',
          text = paste("An error occurred while processing the data:", e$message),
          type = "error"
        )
        return(NULL)
      }, finally = {
        shinyjs::delay(ms = 2000, {
          shinybusy::remove_modal_spinner()
        })
      })

      result
    })

    observe({
     # input$config
      req(Result(),input$parents, Genotype_names())
      # Locking  parents selection
      values$locked_parents <- input$parents
      #print(Result())
    })

    # Get colnames of Mapfile
    map_file_col <- reactive({
      req(Result()$mapfile)
      colnames(Result()$mapfile)
    })



   # Obtain the order of labels
    color_code_checker_res <- reactive({
      req(Result() , values$locked_parents)
      color_code_checker(Result()$proc_kasp_f ,
                         par_1 = values$locked_parents[1],
                         par_2 = values$locked_parents[2])
    })

    # Output the generated color code as a formatted string in the UI
    output$color_format <- renderText({
      req(color_code_checker_res())
      paste(color_code_checker_res(), collapse = ' --> ')

    })


    # Observer for updating map-related inputs
    observe({
      req(map_file_col())

      # Update SNP ID selection
      updateSelectInput(session,
                        inputId = "snp_ids",
                        choices = map_file_col(),
                        selected = safe_grep_match(pattern = 'id',
                                                   choices = map_file_col())

      )

      # Update Chromosome selection
      updateSelectInput(session,
                        inputId = "chr",
                        choices = map_file_col(),
                        selected = safe_grep_match(pattern = 'chr',
                                                   choices = map_file_col())
      )

      # Update Position selection
      updateSelectInput(session,
                        inputId = "chr_pos",
                        choices = map_file_col(),
                        selected = safe_grep_match(pattern = 'pos',
                                                   choices = map_file_col())
      )

      req(color_code_checker_res())
      # Update Position selection
      updateSelectInput(session,
                        inputId = "col_labels",
                        choices = color_code_checker_res() |> unname(),
                        selected = color_code_checker_res() |> unname()
      )

    })

    # Get unique chromosomes present
    uniqchr <- reactive({
      req(input$chr, Result()$mapfile)
      unique(Result()$mapfile[[input$chr]])
    })

   # Update selector input for window or individual chromosome view
    observe({
      updateSelectInput(session, inputId = "chr_vw", choices = uniqchr())
    })


    # Get results for no annotation heatmap
    no_annotate <- reactive({
      req(
        length(values$locked_parents) >= 2, input$chr, input$chr_pos, input$group_sz,
        input$legend_title, input$text_size, input$snp_ids, Result(),
        all(values$locked_parents %in% rownames(Result()$proc_kasp_f))
      )

      col_mapping <- input$col_mapping
      col_labels  <- input$col_labels

      # Wait until both inputs are non null and same length or safely fallback
      if (!is.null(col_mapping) && !is.null(col_labels) && length(col_mapping) >= length(col_labels)) {
        tryCatch({
          cross_qc_heatmap(
            col_mapping = col_mapping[seq_len(length(col_labels))],
            col_labels  = col_labels,
            x = Result()$proc_kasp_f,
            map_file = Result()$mapfile,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents,
            group_sz = if (input$group_sz == 0) nrow(Result()$proc_kasp_f) - 2 else input$group_sz,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
        }, error = function(e) {
          return(NULL)
          shinyWidgets::show_toast(title = '',type = 'error',
                                   text = e$message)
        })


      } else {
        tryCatch({
          cross_qc_heatmap(
            x = Result()$proc_kasp_f,
            map_file = Result()$mapfile,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents,
            group_sz = if (input$group_sz == 0) nrow(Result()$proc_kasp_f) - 2 else input$group_sz,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
        }, error = function(e) {
          shinyWidgets::show_toast(title = '',type = 'error',text = e$message)
          return(NULL)
        })

      }
    })


    # Find the number of batches heatmap has been drawn for
    observe({
      if (!is.null(no_annotate())) {
        updateSelectInput(session, inputId = "batch_no", choices = names(no_annotate()))
      }
    })

    # Update map dynamically.
    heatmap_1 <- reactive({
      req(no_annotate(), input$batch_no)
      if (input$batch_no %in% names(no_annotate())) {
        no_annotate()[[input$batch_no]]
      } else {
        NULL
      }
    })

    # Print plot
    output$heatmap <- renderPlot({
      tryCatch({
        plot_obj <- heatmap_1()
        req(plot_obj)

        if (inherits(plot_obj, "ggplot")) {
          print(plot_obj)
       shinyjs::delay(ms = 1500,{
         shinyWidgets::show_toast(
           title = '',
           text = 'Heatmap Generated Successfully',
           type = "success"
         )
       })

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
        req(no_annotate())  # Ensure plots exist

        tryCatch(
          {
            # Start PDF (onefile=TRUE ensures multi-page)
            grDevices::pdf(file,
                           width = input$width,
                           height = input$height,
                           onefile = TRUE)

            # Print ALL plots regardless of type
            for (plot_obj in no_annotate()) {
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





    ### Annotation Section  ###

    # Store our data and counter
    values <- reactiveValues(
      count = 0, # number of traits
      trait_data = NULL, # Final trait position list
      removed_traits = c() # Track which trait IDs have been removed
    )

    # Function to create one trait input row
    makeTraitInput <- function(id) {
      div(
        id = ns(paste0("trait_", id)),
        style = "border: 1px solid #ddd; padding: 10px; margin: 5px;",
        fluidRow(
          column(3, textInput(ns(paste0("name_", id)), "Name:", value = paste0("loc", id))),
          column(3, numericInput(ns(paste0("chr_", id)), "Chr:", value = 1, min = 1)),
          column(3, numericInput(ns(paste0("start_", id)), "Start:", value = 1000000)),
          column(3, numericInput(ns(paste0("end_", id)), "End:", value = 1400000))
        ),

        # Remove button: only show for traits after the first one
        if (id > 1) {
          actionButton(ns(paste0("remove_", id)), "Remove", class = "btn-danger btn-sm")
        }
      )
    }

    # Create trait is clicked
    observeEvent(input$newBtn, {
      values$count <- 1
      values$removed_traits <- c() # Reset removed traits list

      showModal(modalDialog(
        title = "Define Trait Positions",
        size = "l",

        # Container for dynamic trait inputs
        div(
          id = ns("inputs"),
          makeTraitInput(1)
        ),

        # Footer with action buttons
        footer = tagList(
          actionButton(ns("addTrait"), "Add Trait", class = "btn-success", icon = icon("plus")),
          modalButton("Cancel"),
          actionButton(ns("submit"), "Submit", class = "btn-primary", icon = icon("check"))
        ),
        easyClose = FALSE,
        fade = TRUE
      ))
    })

    # Add new trait when Add Trait is clicked
    observeEvent(input$addTrait, {
      values$count <- values$count + 1
      insertUI(
        selector = paste0("#", ns("inputs")),
        where = "beforeEnd",
        ui = makeTraitInput(values$count)
      )
    })

    # Handle remove buttons dynamically
    observe({
      # Check for remove button clicks for each trait
      for (i in 2:values$count) {
        local({
          my_i <- i # Capture the current i value
          observeEvent(input[[paste0("remove_", my_i)]],
                       {
                         # Remove this trait's UI
                         removeUI(selector = paste0("#", ns(paste0("trait_", my_i))))
                         # Track that this trait ID has been removed
                         values$removed_traits <- c(values$removed_traits, my_i)
                       },
                       ignoreInit = TRUE
          )
        })
      }
    })


    # if batch changes render annotated heatmap as null.
    observeEvent(input$batch, {
      values$trait_data <- NULL
    })


    observeEvent(input$submit, {
      trait_list <- list()

      # Loop through all trait inputs and collect values
      for (i in 1:values$count) {
        # Skip if this trait was removed
        if (i %in% values$removed_traits) {
          next
        }

        name <- input[[paste0("name_", i)]]
        chr <- input[[paste0("chr_", i)]]
        start <- input[[paste0("start_", i)]]
        end <- input[[paste0("end_", i)]]

        # Only add if all values exist and are valid
        if (!is.null(name) && !is.null(chr) && !is.null(start) && !is.null(end) &&
            !is.na(name) && !is.na(chr) && !is.na(start) && !is.na(end) &&
            nchar(trimws(name)) > 0) {
          # Create a named numeric vector (same structure as hardcoded)
          trait_vector <- c(
            chr = as.numeric(chr),
            start = as.numeric(start),
            end = as.numeric(end)
          )

          # Add to list with the name as key
          trait_list[[trimws(name)]] <- trait_vector
        }
      }

      # Store the result and close modal
      values$trait_data <- trait_list
      removeModal()
    })


    # Process trait data with genotype and map files
    neat_result <- reactive({
      req(values$trait_data, Result()$mapfile, Result()$proc_kasp_f)

      tryCatch({
        combi_func(
          list_loc = values$trait_data,
          map_doc = Result()$mapfile,
          genofile = Result()$proc_kasp_f
        )
      }, error = function(e) {
          shinyWidgets::show_alert(
            title = "Error",
            text = paste("Error processing trait data:", e$message),
            type = "error"
          )
        NULL
      })
    })



    # Generate annotated heatmap
    annotated <- reactive({
      req(
        length(values$locked_parents) >= 2, input$chr, input$chr_pos,
        input$text_scale_fct, input$legend_title,
        input$text_size, input$snp_ids,
        values$trait_data, neat_result()
      )

      # hold col mapping and labels
      col_mapping <- input$col_mapping
      col_labels  <- input$col_labels

      # Wait until both inputs are non null and same length  or safely fallback
      if (!is.null(col_mapping) && !is.null(col_labels) && length(col_mapping) >= length(col_labels)) {
        tryCatch({
          cross_qc_heatmap(
            col_mapping = col_mapping[seq_len(length(col_labels))],
            col_labels  = col_labels,
            x = neat_result()$chr_a,
            map_file = neat_result()$chr_map,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents[1:2],
            trait_pos = values$trait_data,
            text_scale_fct = input$text_scale_fct,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
        }, error = function(e) {
          shinyWidgets::show_alert(title = 'Error',
                                   text = e$message,
                                   type = 'error')
          return(NULL)
        })


      } else {
        tryCatch({
          cross_qc_heatmap(
            x = neat_result()$chr_a,
            map_file = neat_result()$chr_map,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents[1:2],
            trait_pos = values$trait_data,
            text_scale_fct = input$text_scale_fct,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
        }, error = function(e) {
          shinyWidgets::show_alert(title = 'Error',
                                   text = e$message,
                                   type = 'error')
          return(NULL)
        })

      }

    })

    # Generate annotated heatmap with second option.
    annotated_2 <- reactive({
      req(
        length(values$locked_parents) >= 2, input$chr, input$chr_pos,
        input$text_scale_fct, input$legend_title,
        input$text_size, input$snp_ids,
        values$trait_data, neat_result()
      )


      col_mapping <- input$col_mapping
      col_labels  <- input$col_labels

      # Wait until both inputs are non null and same length or safely fallback
      if (!is.null(col_mapping) && !is.null(col_labels) && length(col_mapping) >= length(col_labels)) {
        tryCatch({
          cross_qc_heatmap2(
            col_mapping = col_mapping[seq_len(length(col_labels))],
            col_labels  = col_labels,
            x = neat_result()$chr_a,
            map_file = neat_result()$chr_map,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents[1:2],
            trait_pos = values$trait_data,
            text_scale_fct = input$text_scale_fct,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
        }, error = function(e) {
          shinyWidgets::show_alert(title = 'Error',
                                   text = e$message,
                                   type = 'error')
          return(NULL)

        })


      } else {
        tryCatch({
          cross_qc_heatmap2(
            x = neat_result()$chr_a,
            map_file = neat_result()$chr_map,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents[1:2],
            trait_pos = values$trait_data,
            text_scale_fct = input$text_scale_fct,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
        }, error = function(e) {
         shinyWidgets::show_alert(title = 'Error',
                                  text = e$message,
                                  type = 'error')
          return(NULL)
        })

      }
    })




    # Display annotated.
    output$ant_heatmap <- renderPlot({
      if(input$options == '1'){
        req(annotated())
        plot_result <- annotated()

        # Show toast after successful plot generation
        shinyWidgets::show_toast(
          title = '',
          text = 'Box Annotated Heatmap Generated Successfully',
          type = 'success'
        )

        print(plot_result)

      } else if(input$options == '2'){
        req(annotated_2())
        plot_result <- annotated_2()

        # Show toast after successful plot generation
        shinyWidgets::show_toast(
          title = '',
          text = 'Line Annotated Heatmap Generated Successfully',
          type = 'success'
        )

        print(plot_result)
      }
    })

    # Download batch plot.
    output$download_plot2 <- downloadHandler(
      filename = function() {
        req(input$file_name2)
        paste0(input$file_name2, ".pdf")  # Single PDF containing ALL plots
      },
      content = function(file) {
        req(no_annotate())  # Ensure plots exist

        tryCatch(
          {
            # Start PDF (onefile=TRUE ensures multi-page)
            grDevices::pdf(file,
                           width = input$width2,
                           height = input$height2,
                           onefile = TRUE)

            # Print ALL plots regardless of type
            for (plot_obj in list(annotated(),annotated_2())) {
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



    # Window size
    windowsize <- reactive({
      req(input$chr_vw, Result())
      window_size_func(
        data = Result()$proc_kasp_f,
        mapfile = Result()$mapfile,
        chr = input$chr_vw
      )
    })

    #
    wind_result <- reactiveVal()
    observe({
      req(
        windowsize(),
        length(values$locked_parents) >= 2, input$chr, input$chr_pos, # input$group_sz,
        input$legend_title, input$text_size, input$snp_ids, Result()
      )

      col_mapping <- input$col_mapping
      col_labels  <- input$col_labels

      # Wait until both inputs are non null and same length  or safely fallback
      if (!is.null(col_mapping) && !is.null(col_labels) && length(col_mapping) >= length(col_labels)) {
        tryCatch({
        result <-  cross_qc_heatmap(
            col_mapping = col_mapping[seq_len(length(col_labels))],
            col_labels  = col_labels,
            x = windowsize(),
            map_file = Result()$mapfile,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents,
            # trait_pos = values$trait_data,
            # text_scale_fct = input$text_scale_fct,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
        wind_result(result)
        }, error = function(e) {
          shinyWidgets::show_alert(title = 'Error',
                                   text = paste("Heatmap generation error: ", e$message),
                                   type = 'error')
          return(NULL)
        })


      } else {
        tryCatch({
         result <-  cross_qc_heatmap(
            x = windowsize(),
            map_file = Result()$mapfile,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = values$locked_parents,
            # trait_pos = values$trait_data,
            # text_scale_fct = input$text_scale_fct,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            panel_fill = input$panel_fill,
            panel_col = input$panel_col,
            pdf = FALSE
          )
         wind_result(result)
        }, error = function(e) {
         # message("Heatmap generation error: ", e$message)
          return(NULL)
        })

      }

    })

    # Render plot for window size
    output$wind_heatmap <- renderPlot({
      req(wind_result())
      print(wind_result())
    })


    #------------------------------------------------------
    # Information to the user
    #------------------------------------------------------
    # parent missing
    output$par_missing_tbl <- DT::renderDT({
      req(Result())
      DT::datatable(par_missing_dat(Result()$par_missing_dat), options = list(scrollX = TRUE))
    })


    # Genotype error
    output$geno_error_tbl <- DT::renderDT({
      req(Result())
      DT::datatable(Result()$genotype_error, options = list(scrollX = TRUE))
    })

    # # parent hetero
    output$par_het_tbl <- DT::renderDT({
      req(Result())
      DT::datatable(Result()$parent_het, options = list(scrollX = TRUE))
    })

  })
}


