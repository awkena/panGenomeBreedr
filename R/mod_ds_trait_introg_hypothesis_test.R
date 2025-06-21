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
mod_ds_trait_introg_hypothesis_test_ui <- function(id) {
  ns <- NS(id)
  tagList(
    navset_card_underline(
      nav_panel(
        title = "Trait Introgression Hypothesis Testing",
        icon = icon("flask"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            width = 350,
            position = "left",
            bslib::card(
              height = "100%",
              bslib::card_header("Upload & Configure Data"), # class = "bg-primary"),
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
                ),
                radioButtons(
                  inputId = ns("choice"),
                  label = "Do you have a map file?",
                  choices = c("Yes" = "yes", "No, generate one for me" = "no"),
                  selected = "yes"
                ),

                # Dynamic rendering based on choice
                conditionalPanel(
                  condition = paste0('input["', ns("choice"), '"] == "yes"'),
                  tagList(
                    fileInput(
                      inputId = ns("mapfile"),
                      label = "Upload Map file",
                      accept = ".csv",
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
                ),

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
              ),
              tagList(
                bslib::input_switch(
                  id = ns("apply_par_poly"),
                  label = "Remove Polymorphic Parents",
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
              ),
              card_footer(
                actionButton(
                  inputId = ns("config"),
                  label = "Submit",
                  icon = icon("check"),
                  width = "100%",
                  class = "btn-success"
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
                  max_height = "850px",
                  height = "100%",
                  bslib::card_header("Heatmap Configuration", class = "bg-success"),
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
                  max_height = "850px",
                  height = "100%",
                  bslib::card_header("Heatmap Visualization Controls", class = "bg-info"),
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
                          label = "Choose Heatmap Colors (6 Genotypes)",
                          choices = grDevices::colors(),
                          multiple = TRUE,
                          width = "100%"
                        ),
                        selectInput(
                          inputId = ns("col_labels"),
                          label = "Genotype Labels (comma-separated)",
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
                          )
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
                          actionButton(
                            inputId = ns("newBtn"),
                            label = "Create Trait Positions",
                            icon = icon("plus-circle"),
                            class = "btn-primary mb-3",
                            width = '30%'
                          ),
                          radioButtons(
                            inputId = ns("options"),
                            label = "View Heatmap Options",
                            choices = c(
                              "Options 1" = "1",
                              "Options 2" = "2"
                            ),
                            selected = "1",
                            inline = TRUE,
                            width = "500px"
                          ),
                          plotOutput(
                            outputId = ns("ant_heatmap"),
                            width = "100%",
                            height = "600px"
                          )
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
            title = "📌 Important Note",
            tagList(
              p("Marker names must follow a structured format to be parsed into the map file."),
              tags$ul(
                tags$li("A common prefix before each marker (e.g., 'S')."),
                tags$li("Chromosome number immediately after the prefix (e.g., '1')."),
                tags$li("A separator character (e.g., '_')."),
                tags$li("Position number after the separator (e.g., '101').")
              ),
              p("Example: 'S1_101' — where 'S' is the prefix, '1' is the chromosome number, and '101' is the position.")
            ),
            easyClose = FALSE,
            footer = modalButton("Got it!")
          )
        )
      }
    })


    # Read the csv file uploadeed.
    data <- reactive({
      req(input$data_id)
      read.csv(file = input$data_id$datapath) |> as.data.frame()
    })

    # Get colnames of data and populate.
    data_colnames <- reactive({
      req(data())
      Get_dt_coln(data())
    })

    # Populate batch and genotype columns with data colnames
    observe({
      updateSelectInput(session,
        inputId = "batch_col",
        choices = data_colnames(),
        selected = grep("batch",
          x = data_colnames(),
          ignore.case = TRUE,
          value = TRUE
        )[1]
      )

      updateSelectInput(session,
        inputId = "genotype_col",
        choices = data_colnames(),
        selected = grep("genotype",
          x = data_colnames(),
          ignore.case = TRUE,
          value = TRUE
        )[1]
      )
    })




    # Get unique batch from it.
    uniq_batch <- reactive({
      req(data(), input$batch_col)
      data()[[input$batch_col]] |> unique()
    })

    # Populate batch widget
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
      req(input$batch, data())
      Genotypes_user(data = data(), Batch = input$batch )
    })

    observe({
      req(Genotype_names())
      updateSelectInput(session, inputId = "dp", choices = Genotype_names())
      updateSelectInput(session, inputId = "rp", choices = Genotype_names(), selected = Genotype_names()[3])
    })

    # Read map file if user has.
    map_file <- reactive({
    req(input$mapfile)
    read_mapfile(file = input$mapfile$datapath)

    })

    # Populate field for snp id column based on mapfile
    observe({
      req(map_file())
        updateSelectInput(session,
                          inputId = "snp_id",
                          choices = colnames(map_file()),
                          selected = grep("id",
                                          x = colnames(map_file()), ignore.case = TRUE,
                                          value = TRUE)[1]
        )
    })


    # store the parents selection temporary
    values <- reactiveValues(
      locked_parents = NULL
    )


    Result <- reactiveVal(NULL) # empty reactiveval object

    # Generate list of mapfile, proccess data etc when submit is clicked
    observeEvent(input$config, {
      req(
        data(), input$batch, input$batch_col,
        input$data_type, input$allele_sep,
        input$rp, input$dp, input$choice,
        input$genotype_col
      )

      # Process the data
      result <- proc_nd_map_func(
        data = data(),
        Batch = input$batch,
        batch_col = input$batch_col,
        marker_sep = if(input$choice == "no") input$sep_marker else NULL,
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
        Prefix = if (input$choice == "no") isolate(input$prefix_marker) else NULL,
        geno_vec = Genotype_names(),
        feedback = input$choice,
        na_code = NA,
        data_col = data_colnames(),
        mapfile_path = map_file()
      )

      # Store result
      Result(result)

      # Update parents selection
      req(Genotype_names())

      updateSelectInput(session,
                        inputId = "parents",
                        choices = Genotype_names(),
                        selected = Genotype_names()[c(1, 3)]
      )

      # Locking  parents selection
      values$locked_parents <- Genotype_names()[c(1, 3)]
    })



    # Get colnames of Mapfile
    map_file_col <- reactive({
      req(Result()$mapfile)
      colnames(Result()$mapfile)
    })



    # Update the verbatim print out put.

    color_code_checker_res <- reactive({
      req(Result())
      color_code_checker(Result()$proc_kasp_f)
    })


    output$color_format <- renderText({
      req(color_code_checker_res())
      paste(names(color_code_checker_res()$color_code_in), collapse = ' --> ')
    })


    # Observer for updating map-related inputs
    observe({
      req(Result()$mapfile, map_file_col())

      # Update SNP ID selection
      updateSelectInput(session,
                        inputId = "snp_ids",
                        choices = map_file_col(),
                        selected = grep(
                          pattern = "id", x = map_file_col(),
                          ignore.case = TRUE, value = TRUE
                        )[1]
      )

      # Update Chromosome selection
      updateSelectInput(session,
                        inputId = "chr",
                        choices = map_file_col(),
                        selected = grep(
                          pattern = "chr", x = map_file_col(),
                          ignore.case = TRUE, value = TRUE
                        )[1]
      )

      # Update Position selection
      updateSelectInput(session,
                        inputId = "chr_pos",
                        choices = map_file_col(),
                        selected = grep(
                          pattern = "pos", x = map_file_col(),
                          ignore.case = TRUE, value = TRUE
                        )[1]
      )

      req(color_code_checker_res())
      # Update Position selection
      updateSelectInput(session,
                        inputId = "col_labels",
                        choices = color_code_checker_res()$just_num,
                        selected = color_code_checker_res()$just_num
      )

    })

    uniqchr <- reactive({
      req(input$chr, Result()$mapfile)
      unique(Result()$mapfile[[input$chr]])
    })


    observe({
      updateSelectInput(session, inputId = "chr_vw", choices = uniqchr())
    })







    #Using locked parents in heat map to prevent crashes when batch number changes
    no_annotate <- reactive({
      req(
        values$locked_parents, input$chr, input$chr_pos, input$group_sz,
        input$legend_title, input$text_size, input$snp_ids, Result()
      )

      # Validate locked parents selection
      if (length(values$locked_parents) < 2) {
        return(NULL)
      }

      cross_qc_heatmap(
        x = Result()$proc_kasp_f,
        map_file = Result()$mapfile,
        snp_ids = input$snp_ids,
        chr = input$chr,
        chr_pos = input$chr_pos,
        parents = values$locked_parents, # Use LOCKED parents, not input$parents
        group_sz = if (input$group_sz == 0) nrow(Result()$proc_kasp_f) - 2 else input$group_sz,
        legend_title = input$legend_title,
        text_size = input$text_size,
        alpha = input$alpha,
        panel_fill = input$panel_fill,
        panel_col = input$panel_col,
        pdf = FALSE
      )
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

    # Render default heatmap
    output$heatmap <- renderPlot({
      req(heatmap_1())
      print(heatmap_1())
    })





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

      observeEvent(input$submit, {
      # shinybusy::show_modal_spinner(
      #   spin = "fading-circle",
      #   color = "#0dc5c1",
      #   text = "Generating Annotated Heatmap... Please wait."
      # )
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
        if (exists("shinyWidgets") && "show_alert" %in% names(getNamespace("shinyWidgets"))) {
          shinyWidgets::show_alert(
            title = "Error",
            text = paste("Error processing trait data:", e$message),
            type = "error"
          )
        }
        NULL
      })
    })



    # Generate annotated heatmap
    annotated <- reactive({
      req(
        input$parents, input$chr, input$chr_pos,
        input$text_scale_fct, input$legend_title,
        input$text_size, input$snp_ids,
        values$trait_data, neat_result()
      )

      # Validate parents
      if (length(input$parents) < 2) {
        return(NULL)
      }

      tryCatch({
        cross_qc_heatmap(
          x = neat_result()$chr_a,
          map_file = neat_result()$chr_map,
          snp_ids = input$snp_ids,
          chr = input$chr,
          chr_pos = input$chr_pos,
          parents = input$parents[1:2],
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
        if (exists("shinyWidgets") && "show_toast" %in% names(getNamespace("shinyWidgets"))) {
          shinyWidgets::show_toast(
            title = "",
            text = paste("Error in cross_qc_heatmap:", e$message),
            type = "error"
          )
        }
        NULL
      })
    })


    # Generate annotated heatmap
    annotated_2 <- reactive({
      req(
        input$parents, input$chr, input$chr_pos,
        input$text_scale_fct, input$legend_title,
        input$text_size, input$snp_ids,
        values$trait_data, neat_result()
      )

      # Validate parents
      if (length(input$parents) < 2) {
        return(NULL)
      }

      tryCatch({
        cross_qc_heatmap2(
          x = neat_result()$chr_a,
          map_file = neat_result()$chr_map,
          snp_ids = input$snp_ids,
          chr = input$chr,
          chr_pos = input$chr_pos,
          parents = input$parents[1:2],
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
        if (exists("shinyWidgets") && "show_toast" %in% names(getNamespace("shinyWidgets"))) {
          shinyWidgets::show_toast(
            title = "",
            text = paste("Error in cross_qc_heatmap:", e$message),
            type = "error"
          )
        }
        NULL
      })
    })




    # Display annotated.
    observe({
      if(input$options == '1'){
        output$ant_heatmap <- renderPlot({
          req(annotated())
          print(annotated())
        })
      } else if(input$options == '2'){
        output$ant_heatmap <- renderPlot({
          req(annotated_2())
          print(annotated_2())
        })
      }
    })

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

    wind_result <- eventReactive(input$chr_vw, {
      req(
        windowsize(),
        input$parents, input$chr, input$chr_pos, # input$group_sz,
        input$legend_title, input$text_size, input$snp_ids, Result()
      )

      # Validate that selected parents still exist in current genotype names
      current_genotypes <- Genotype_names()
      valid_parents <- input$parents[input$parents %in% current_genotypes]

      # Ensure we have at least 2 valid parents
      if (length(valid_parents) < 2) {
        # Show alert
        shinyWidgets::show_alert(
          title = "Error",
          text = "Select a Recurrent & Donor Parent",
          type = "error"
        )

        return(NULL)
      }

      # Use only the first 2 valid parents
      parents_to_use <- valid_parents[1:2]

      cross_qc_heatmap(
        x = windowsize(),
        map_file = Result()$mapfile,
        snp_ids = input$snp_ids,
        chr = input$chr,
        chr_pos = input$chr_pos,
        parents = parents_to_use,
        # trait_pos = values$trait_data,
        # text_scale_fct = input$text_scale_fct,
        legend_title = input$legend_title,
        text_size = input$text_size,
        alpha = input$alpha,
        panel_fill = input$panel_fill,
        panel_col = input$panel_col,
        pdf = FALSE
      )
    })

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
      DT::datatable(Result()$par_missing_dat, options = list(scrollX = TRUE))
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


