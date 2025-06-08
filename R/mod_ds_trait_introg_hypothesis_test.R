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
              bslib::card_header("Upload & Configure Data", class = "bg-primary"),
              bslib::card_body(
                fileInput(
                  inputId = ns("data_id"),
                  label = "Upload Kasp/Agriplex File",
                  multiple = FALSE,
                  accept = ".csv",
                  width = "100%"
                ),
                radioButtons(
                  inputId = ns("choice"),
                  label = "Do you have a map file?",
                  choices = c("Yes" = "yes", "No, generate one for me" = "no"),
                  selected = "no"
                ),
                uiOutput(outputId = ns("map_upld")),
                selectInput(
                  inputId = ns("batch"),
                  label = "Select Focused Batch",
                  choices = NULL,
                  width = "100%"
                ),
                textInput(
                  inputId = ns("sep_marker"),
                  label = "Enter Separator for Marker ID",
                  value = "_",
                  width = "100%"
                ),
                selectInput(
                  inputId = ns("data_type"),
                  label = "Indicate Data Format",
                  choices = c("agriplex", "Kasp"),
                  selected = "agriplex",
                  width = "100%"
                ),
                textInput(
                  inputId = ns("allele_sep"),
                  label = "Enter Allele Separator for Data Type",
                  value = " / ",
                  width = "100%"
                ),
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
          fluidRow(
            # Card 1: Heatmap Configuration
            column(
              width = 4,
              bslib::card(
                max_height = "750px",
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
                    value = 1,
                    min = 1,
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
                max_height = "750px",
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
                      textInput(
                        inputId = ns("col_labels"),
                        label = "Genotype Labels (comma-separated)",
                        value = NULL,
                        width = "100%"
                      ),
                      selectInput(
                        inputId = ns("panel_fill"),
                        label = "Panel Background Fill Color",
                        choices = grDevices::colors(),
                        selected = "grey80",
                        width = "100%"
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
                      title = "Basic Heatmaps",
                      icon = icon("th"),
                      fluidRow(
                        column(
                          width = 4,
                          selectInput(
                            inputId = ns("batch_no"),
                            label = "Select Batch for Display",
                            choices = NULL,
                            width = "100%"
                          )
                        )
                      ),
                      plotOutput(
                        outputId = ns("heatmap"),
                        width = "100%",
                        height = "600px"
                      )
                    ),
                    bslib::nav_panel(
                      title = "Annotated Heatmaps",
                      icon = icon("tags"),
                      fluidRow(
                        column(
                          width = 12,
                          fluidRow(
                            column(
                              width = 4,
                              numericInput(
                                inputId = ns("mk_plot"),
                                label = "Number of Markers to Plot",
                                value = 20,
                                min = 1
                              )
                            ),
                            column(
                              width = 4,
                              numericInput(
                                inputId = ns("text_scale_fct"),
                                label = "Text Scaling Size",
                                value = 0.18,
                                min = 0.1,
                                max = 1
                              )
                            )
                          ),
                          actionButton(
                            inputId = ns("newBtn"),
                            label = "Create Trait Positions",
                            icon = icon("plus-circle"),
                            class = "btn-primary mb-3"
                          )
                        )
                      ),
                      plotOutput(
                        outputId = ns("heatmap_ann"),
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

    # Dynamic ui based on choice of individual.
    observe({
      req(input$choice)
      if (input$choice == "yes") {
        output$map_upld <- renderUI({
          fileInput(
            inputId = ns("mapfile"),
            label = "Upload Map file",
            accept = ".csv",
            width = "100%"
          )
        })
      } else if (input$choice == "no") {
        output$map_upld <- renderUI({
          NULL
        })
      }
    })

    # Read the csv file uploadeed.
    data <- reactive({
      req(input$data_id)
      read.csv(file = input$data_id$datapath) |> as.data.frame()
    })

    # Get unique batch from it.
    uniq_batch <- reactive({
      req(data())
      batch <- grep(pattern = "batch", x = colnames(data()), ignore.case = TRUE, value = TRUE)
      data()[[batch]] |> unique()
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
      Genotypes_user(data = data(), Batch = input$batch)
    })

    observe({
      req(Genotype_names())
      updateSelectInput(session, inputId = "dp", choices = Genotype_names())
      updateSelectInput(session, inputId = "rp", choices = Genotype_names())
    })

    # Read map file if user has.
    map_file <- reactive({
      if (is.null(input$mapfile)) {
        NULL
      } else {
        read.csv(file = input$mapfile$datapath)
      }
    })

    #-- Allow users to submit.
    Result <- reactiveVal(NULL) # empty reactive value to store result.

    observeEvent(input$config, {
      req(
        data(), input$batch, input$sep_marker, Genotype_names(),
        input$data_type, input$allele_sep, input$rp, input$dp, input$choice
      )
      # Cleaning and numeric coding
      result <- marker_file(
        data = data(),
        Batch = input$batch,
        sep = input$sep_marker,
        sep_2 = check_sep(input$allele_sep),
        data_type = input$data_type,
        rp = input$rp,
        dp = input$dp,
        geno_vec = Genotype_names(),
        feedback = input$choice,
        mapfile_path = if (is.null(map_file())) NULL else map_file()
      )

      # Store result in reactive value
      Result(result)

      # Reset the parents selection to avoid stale selections
      updateSelectInput(session,
        inputId = "parents",
        choices = Genotype_names(),
        selected = NULL
      )

      # Then set default selection after a brief delay
      Sys.sleep(0.1)
      updateSelectInput(session,
        inputId = "parents",
        choices = Genotype_names(),
        selected = Genotype_names()[c(1, min(3, length(Genotype_names())))]
      )
    })

    # Get colnames of Mapfile
    map_file_col <- reactive({
      req(Result()$mapfile)
      colnames(Result()$mapfile)
    })

    # observer for updating map-related inputs
    observeEvent(Result(), {
      req(Result()$mapfile, map_file_col(), Genotype_names())

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
    })


    # Draw the heat map
    # Add validation for parents selection
    no_annotate <- reactive({
      req(
        input$parents, input$chr, input$chr_pos, input$group_sz,
        input$legend_title, input$text_size, input$snp_ids, Result()
      )

      # Validate that selected parents still exist in current genotype names
      current_genotypes <- Genotype_names()
      valid_parents <- input$parents[input$parents %in% current_genotypes]

      # Ensure we have at least 2 valid parents
      if (length(valid_parents) < 2) {
        return(NULL)
      }

      # Use only the first 2 valid parents
      parents_to_use <- valid_parents[1:2]

      cross_qc_heatmap(
        x = Result()$proc_kasp_f,
        map_file = Result()$mapfile,
        snp_ids = input$snp_ids,
        chr = input$chr,
        chr_pos = input$chr_pos,
        parents = parents_to_use, # Use validated parents
        group_sz = input$group_sz,
        legend_title = input$legend_title,
        text_size = input$text_size,
        alpha = input$alpha,
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

    # Render plot
    output$heatmap <- renderPlot({
      if (!is.null(heatmap_1())) {
        print(heatmap_1())
      }
    })


    ### ---- Annotation aspect.

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

    # Fixed version of your submit handler to ensure proper data structure
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


    # Updated annotated reactive with better error handling
    annotated <- reactive({
      req(
        input$parents, input$chr, input$chr_pos, input$mk_plot, input$text_scale_fct,
        input$legend_title, input$text_size, input$snp_ids, Result(), values$trait_data
      )

      # Validate that selected parents still exist in current genotype names
      current_genotypes <- Genotype_names()
      valid_parents <- input$parents[input$parents %in% current_genotypes]

      # Ensure we have at least 2 valid parents
      if (length(valid_parents) < 2) {
        showNotification("Less than 2 valid parents selected", type = "warning")
        return(NULL)
      }

      # Use only the first 2 valid parents
      parents_to_use <- valid_parents[1:2]

      # Validate trait_data structure
      if (is.null(values$trait_data) || length(values$trait_data) == 0) {
        showNotification("No trait data available", type = "error")
        return(NULL)
      }

      # Check that trait_data has the expected structure
      for (trait_name in names(values$trait_data)) {
        trait_info <- values$trait_data[[trait_name]]
        if (!all(c("chr", "start", "end") %in% names(trait_info))) {
          showNotification(paste("Invalid trait data structure for:", trait_name), type = "error")
          return(NULL)
        }
        if (any(is.na(trait_info))) {
          showNotification(paste("NA values in trait data for:", trait_name), type = "error")
          return(NULL)
        }
      }

      # Subset the mapfile
      uniq_chr <- sapply(values$trait_data, function(x) x["chr"]) |> unique()
      if (length(uniq_chr) > 1) {
        showNotification("All chromosomes selected must be the same", type = "error")
        return(NULL)
      }

      if (length(uniq_chr) == 0 || is.na(uniq_chr)) {
        showNotification("No valid chromosome found in trait data", type = "error")
        return(NULL)
      }

      # Get chromosome subset
      chr_logical <- Result()$mapfile["chr"] == uniq_chr
      if (sum(chr_logical, na.rm = TRUE) == 0) {
        showNotification("No markers found for selected chromosome", type = "error")
        return(NULL)
      }

      # Get the subset data - use the same data for both x and the chr_a calculation
      chr_a <- Result()$proc_kasp_f[, chr_logical, drop = FALSE]
      max_markers <- min(input$mk_plot, ncol(chr_a))

      if (max_markers <= 0) {
        showNotification("No markers available for plotting", type = "error")
        return(NULL)
      }

      chr_a <- chr_a[, seq_len(max_markers), drop = FALSE]
      chr_map <- parse_marker_ns(colnames(chr_a))

      # Final validation
      if (is.null(chr_map) || nrow(chr_map) == 0) {
        showNotification("Failed to parse marker names", type = "error")
        return(NULL)
      }

      tryCatch(
        {
          cross_qc_heatmap(
            x = chr_a,
            map_file = chr_map,
            snp_ids = input$snp_ids,
            chr = input$chr,
            chr_pos = input$chr_pos,
            parents = parents_to_use,
            trait_pos = values$trait_data,
            text_scale_fct = input$text_scale_fct,
            legend_title = input$legend_title,
            text_size = input$text_size,
            alpha = input$alpha,
            pdf = FALSE
          )
        },
        error = function(e) {
          showNotification(paste("Error in cross_qc_heatmap:", e$message), type = "error")
          return(NULL)
        }
      )
    })


    # Display annotated.
    output$heatmap_ann <- renderPlot({
      req(annotated())
      print(annotated())
    })


    #------------------------------------------------------
    # Information to the user
    #------------------------------------------------------
    # parent missing
    output$par_missing_tbl <- renderDT({
      req(Result())
      DT::datatable(Result()$par_missing_dat, options = list(scrollX = TRUE))
    })

    # Genotype error
    output$geno_error_tbl <- renderDT({
      req(Result())
      DT::datatable(Result()$genotype_error, options = list(scrollX = TRUE))
    })

    # # parent hetero
    output$par_het_tbl <- renderDT({
      req(Result())
      DT::datatable(Result()$parent_het, options = list(scrollX = TRUE))
    })
  })
}
