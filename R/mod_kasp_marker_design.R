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

  # Accordion section title with an icon
  section_title <- function(icon_name, text) {
    div(icon(icon_name, class = "me-2"), tags$span(text, class = "fw-semibold"))
  }

  tagList(
    bslib::layout_sidebar(

      # SIDEBAR ---------------------------------------------------------------

      sidebar = bslib::sidebar(
        id = ns("sidebar"),
        width = 350,
        class = "bg-light",

        bslib::accordion(
          id = "config_accordion",
          open = "files", # Only the first panel is open by default

          # Upload files
          bslib::accordion_panel(
            title = section_title("upload", "Data Acquisition"),
            value = "files",

            # Radio buttons for file type selection
            radioButtons(
              inputId = ns("upload_choice"),
              label = label_help(
                "Variant Data Type",
                "snpEff annotated VCF: the VCF straight from variant calling and annotation. Genotype matrix: a table you have already processed, one row per variant."
              ),
              choices = c("snpEff Annotated VCF", "Genotype Matrix (Processed)"),
              selected = "snpEff Annotated VCF"
            ),
            uiOutput(ns("choice_output")), # dynamic file input based on radio button selection

            # File input for genome reference
            fileInput(
              inputId = ns("genome_file"),
              label = label_help(
                "Reference Genome",
                "FASTA file (.fa, .fasta or .fa.gz) of the reference the variants were called against. A file with only the target chromosome is enough and loads faster."
              ),
              placeholder = ".fa, .fasta, .gz",
              accept = c(".fa", ".fasta", ".gz")
            ),
            example_data_link(ns("load_example")),
            uiOutput(ns("example_badge")),
            tags$details(
              class = "small text-muted mb-3",
              tags$summary("Data format"),
              tags$ul(
                class = "ps-3 mt-2 mb-0",
                tags$li("VCF: a standard VCF annotated with snpEff. Sample genotypes start at column 10."),
                tags$li("Genotype matrix: one row per variant with ID, chromosome, position, reference and alternate allele columns, then one column per sample."),
                tags$li("The chromosome names must match the names in the reference FASTA (e.g. Chr02).")
              )
            )
          ),

          # Map columns
          bslib::accordion_panel(
            title = section_title("project-diagram", "Map Columns"),
            value = "mapping",
            selectizeInput(
              inputId = ns("variant_id_col"),
              label = "Variant ID Column",
              choices = NULL
            ),
            bslib::layout_columns(
              col_widths = c(6, 6),
              selectizeInput(inputId = ns("chrom_col"), label = "Chromosome", choices = NULL),
              selectizeInput(inputId = ns("pos_col"), label = "Position", choices = NULL)
            ),
            bslib::layout_columns(
              col_widths = c(6, 6),
              selectizeInput(inputId = ns("ref_al_col"), label = "Reference Allele", choices = NULL),
              selectizeInput(inputId = ns("alt_al_col"), label = "Alternate Allele", choices = NULL)
            ),
            numericInput(
              inputId = ns("geno_start"),
              label = label_help(
                "Genotypes Start at Column",
                "Number of the first sample column. 10 for a standard VCF."
              ),
              value = 10
            )
          ),

          # Markers to design
          bslib::accordion_panel(
            title = section_title("crosshairs", "Target Variants"),
            value = "markers",
            selectizeInput(
              inputId = ns("chr_ID"),
              label = "Target Chromosome",
              choices = NULL,
              options = list(placeholder = "Choose chromosome...")
            ),
            selectizeInput(
              inputId = ns("marker_ID"),
              label = label_help(
                "Variants to Design",
                "One KASP assay is designed for each variant you pick. All must be on the target chromosome."
              ),
              choices = NULL,
              multiple = TRUE,
              options = list(placeholder = "Select one or more variants...")
            ),
            numericInput(
              inputId = ns("maf"),
              label = label_help(
                "Minor Allele Frequency (0-1)",
                "Other variants near the target whose minor allele frequency is at or above this are marked in the flanking sequence, so the primers can avoid them."
              ),
              value = 0.05, min = 0, max = 1, step = 0.01
            )
          ),

          # Intertek order details
          bslib::accordion_panel(
            title = section_title("file-signature", "Order Details"),
            value = "order",
            tags$p(class = "small text-muted", "These fill the Intertek submission table."),
            textInput(inputId = ns("reg_name"), label = "Region Name", placeholder = "lgs1"),
            textInput(inputId = ns("genome_version"), label = "Genome Version", placeholder = "Sbicolor v5.1"),
            textInput(inputId = ns("trait"), label = "Trait", placeholder = "Stay-green"),
            textInput(inputId = ns("owner"), label = "Owner", placeholder = "Green Evolution")
          )
        ),

        run_button(ns("run_but"), "Design KASP Markers", icon_name = "dna")
      ),

      # PLACEHOLDER -----------------------------------------------------------

      # Shown until markers have been designed
      div(
        id = ns("placeholder_ui"),
        class = "text-center text-muted py-5",
        div(class = "h1 mb-3", icon("dna")),
        tags$p(
          "Upload a variant file and the reference genome in the sidebar, or load the example data, pick the variants, then click",
          tags$b("Design KASP Markers"), "."
        ),
        tags$p(
          class = "small",
          "You'll get the sequence for each assay, ready to send to Intertek, and a plot of each sequence aligned to the reference."
        )
      ),

      # RESULTS ---------------------------------------------------------------

      # Hidden until markers have been designed
      shinyjs::hidden(
        div(
          id = ns("analysis_ui"),
          bslib::navset_card_underline(
            id = ns("results_tabs"),

            # Marker table tab
            bslib::nav_panel(
              title = "Marker Table",
              icon = icon("table"),
              value = "table_tab",
              ui_result_card(
                "KASP assays for Intertek",
                full_screen = TRUE,
                actions = bslib::popover(
                  tags$button(
                    type = "button",
                    class = "btn btn-sm btn-success",
                    icon("download"), " Download"
                  ),
                  textInput(inputId = ns("file_name"), label = "File name", value = "KASP_markers_for_Intertek"),
                  selectInput(
                    inputId = ns("exten"),
                    label = "Format",
                    choices = c("Excel (.xlsx)" = ".xlsx", "CSV (.csv)" = ".csv"),
                    selected = ".xlsx"
                  ),
                  downloadButton(ns("download_table"), "Save", class = "btn-success w-100"),
                  title = "Download table",
                  placement = "bottom"
                ),
                uiOutput(ns("table_summary")),
                DT::DTOutput(ns("kasp_table"))
              )
            ),

            # Alignment plot tab
            bslib::nav_panel(
              title = "Alignment Plot",
              icon = icon("chart-bar"),
              value = "plot_tab",
              ui_result_card(
                "Sequence alignment",
                full_screen = TRUE,
                actions = tagList(
                  div(
                    style = "min-width: 220px;",
                    selectizeInput(
                      inputId = ns("plot_choice"),
                      label = NULL,
                      choices = NULL,
                      width = "100%"
                    )
                  ),
                  downloadButton(
                    outputId = ns("download_plot"),
                    label = "PDF (all)",
                    class = "btn-sm btn-success",
                    icon = icon("download")
                  )
                ),
                tags$p(
                  class = "small text-muted",
                  "The assay sequence for the selected variant aligned to the reference genome around it. Use the PDF button to save the plots for all designed variants."
                ),
                plotOutput(ns("plot"), height = "450px")
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
mod_kasp_marker_design_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # FILE UPLOAD -----------------------------------------------------------

    # VCF or genotype matrix: clear the other input when they switch
    observeEvent(input$upload_choice, {
      if (input$upload_choice == "snpEff Annotated VCF") {
        output$choice_output <- renderUI({
          fileInput(ns("vcf_file"),
            label = "Upload snpEff Annotated VCF",
            accept = c(".vcf",'.gz')
          )
        })
        # Reset existing genotype data input
        shinyjs::reset(ns("gt_df"))
      } else {
        output$choice_output <- renderUI({
          fileInput(ns("gt_df"),
            label = "Upload Genotype Matrix (Processed)",
            accept = c(".xlsx", ".csv")
          )
        })
        # Reset existing VCF file input
        shinyjs::reset(id = ns("vcf_file"))
      }
    })

    # VCF to read: the latest upload, or the bundled example
    vcf_path <- reactiveVal(NULL)
    example_active <- reactiveVal(FALSE)
    example_vcf <- "Sobic.002G302700_SNP_snpeff.vcf"
    example_genome_url <- "https://raw.githubusercontent.com/awkena/panGB/main/Chr02.fa.gz"

    # New VCF upload
    observeEvent(input$vcf_file, {
      vcf_path(input$vcf_file$datapath)
      example_active(FALSE)
    })

    # EXAMPLE DATA ----------------------------------------------------------

    # kasp_marker_design() help page example: a SNP in Sobic.002G302700 on
    # Chr02. The Chr02 genome is not shipped, so we point them to the download
    observeEvent(input$load_example, {
      updateRadioButtons(session, "upload_choice", selected = "snpEff Annotated VCF")
      updateNumericInput(session, "geno_start", value = 10)
      updateNumericInput(session, "maf", value = 0.05)
      updateTextInput(session, "reg_name", value = "ma1")
      updateTextInput(session, "genome_version", value = "Sbv5.1")
      updateTextInput(session, "trait", value = "Example trait")
      updateTextInput(session, "owner", value = "panGenomeBreedr example")
      vcf_path(example_file(example_vcf))
      example_active(TRUE)
    })

    # Clear the example
    observeEvent(input$clear_example, {
      vcf_path(NULL)
      example_active(FALSE)
    })

    # "Using example data" badge
    output$example_badge <- renderUI({
      req(example_active())
      example_data_badge(
        example_vcf,
        clear_id = ns("clear_example"),
        note = tagList(
          tags$strong("One more step: "),
          "the reference genome is not included. ",
          tags$a(href = example_genome_url, target = "_blank", "Download Chr02.fa.gz"),
          " (BTx623 v5.1, chromosome 2), upload it under Reference Genome, then click Design KASP Markers."
        )
      )
    })

    # READ DATA -------------------------------------------------------------

    # Read the VCF
    vcf_data <- reactive({
      req(vcf_path())
      marker.chr_ID(vcf_path())
    })

    # Get column names of VCF file
    vcf_colnames <- reactive({
      req(vcf_path())
      colnames(read_vcf_as_df(vcf_path()))
    })

    # Read the genotype matrix (csv or excel)
    gt_data <- reactive({
      req(input$gt_df)
      tryCatch(
        {
          result <- read_mapfile(filepath = input$gt_df$datapath)

          return(as.data.frame(result))
        },
        error = function(e) {
          shinyWidgets::show_toast(
            title = "Couldn't read the genotype matrix",
            text = e$message,
            type = "error",
            timer = 2000,
            position = "bottom-end"
          )
          return(NULL)
        }
      )
    })

    # Column names of the genotype matrix
    gt_colnames <- reactive({
      req(gt_data())
      # Extract needed column names using safe_grep_match
      safe_grep_match(
        pattern = "variant_id|marker_id|chro|pos|ref|alt",
        choices = colnames(gt_data())
      )
    })

    # Extract unique chromosomes
    unique_chrom <- reactive({
      req(gt_data(), gt_colnames())

      chrm <- safe_grep_match(
        pattern = "chro",
        choices = gt_colnames()
      ) # look for chromosome column
      # Validate chromosome column existence
      if (length(chrm) == 0) {
        shinyWidgets::show_toast(
          title = "",
          text = "No chromosome column found. Name one column chrom or chromosome.",
          type = "error"
        )
        return(NULL)
      }
      # Extract unique chromosomes
      unique(gt_data()[[chrm]])
    })

    # Extract unique marker IDs
    unique_marker_id <- reactive({
      req(gt_data(), gt_colnames())

      mkr <- safe_grep_match(
        pattern = "variant_id|marker_id", # either variant_id or marker_id
        choices = gt_colnames()
      ) # look for variant id column
      # Validate variant id existence
      if (length(mkr) == 0) {
        shinyWidgets::show_toast(
          title = "",
          text = "No variant ID column found. Name one column variant_id or marker_id.",
          type = "error"
        )
        return(NULL)
      }
      # Extract unique variant ids
      unique(gt_data()[[mkr]])
    })

    # COLUMN DROPDOWNS ------------------------------------------------------

    # Fill the dropdowns from the VCF
    observeEvent(vcf_data(), {
      req(vcf_data(), vcf_colnames())

      # Update variable select inputs
      updateSelectizeInput(
        server = TRUE, session, "variant_id_col",
        choices = vcf_colnames(),
        selected = safe_grep_match(
          pattern = "id",
          choices = vcf_colnames()
        )[1]
      )
      # Update chromosome select input
      updateSelectizeInput(
        server = TRUE, session, "chrom_col",
        choices = vcf_colnames(),
        selected = safe_grep_match(
          pattern = "chro",
          choices = vcf_colnames()
        )
      )
      # Update position select input
      updateSelectizeInput(
        server = TRUE, session, "pos_col",
        choices = vcf_colnames(),
        selected = safe_grep_match(
          pattern = "pos",
          choices = vcf_colnames()
        )
      )
      # Update ref allele select inputs
      updateSelectizeInput(
        server = TRUE, session, "ref_al_col",
        choices = vcf_colnames(),
        selected = safe_grep_match(
          pattern = "ref",
          choices = vcf_colnames()
        )
      )
      # Update alt allele select inputs
      updateSelectizeInput(
        server = TRUE, session, "alt_al_col",
        choices = vcf_colnames(),
        selected = safe_grep_match(
          pattern = "alt",
          choices = vcf_colnames()
        )
      )
      # Update marker ID select input allowing multiple selections
      updateSelectizeInput(
        server = TRUE,
        session = session,
        inputId = "marker_ID",
        choices = vcf_data()$vcf_matrix_markerID,
        selected = if (isTRUE(example_active())) "SNP_Chr02_69200443" else NULL
      )
      # Update chromosome ID select input
      updateSelectizeInput(
        server = TRUE,
        session = session,
        inputId = "chr_ID",
        choices = vcf_data()$vcf_matrix_chromID,
        selected = if (isTRUE(example_active())) "Chr02" else NULL
      )
    })

    # Fill the dropdowns from the genotype matrix
    observeEvent(gt_data(), {
      req(unique_marker_id(), unique_chrom(), gt_colnames())
      # Update variable select inputs
      updateSelectizeInput(
        server = TRUE, session, "variant_id_col",
        choices = gt_colnames(),
        selected = safe_grep_match(
          pattern = "id",
          choices = gt_colnames()
        )
      )
      # Update chromosome select input
      updateSelectizeInput(
        server = TRUE, session, "chrom_col",
        choices = gt_colnames(),
        selected = safe_grep_match(
          pattern = "chro",
          choices = gt_colnames()
        )
      )
      # Update position select input
      updateSelectizeInput(
        server = TRUE, session, "pos_col",
        choices = gt_colnames(),
        selected = safe_grep_match(
          pattern = "pos",
          choices = gt_colnames()
        )
      )
      # Update ref allele select input
      updateSelectizeInput(
        server = TRUE, session, "ref_al_col",
        choices = gt_colnames(),
        selected = safe_grep_match(
          pattern = "ref",
          choices = gt_colnames()
        )
      )
      # Update alt allele select input
      updateSelectizeInput(
        server = TRUE,
        session = session,
        inputId = "alt_al_col",
        choices = gt_colnames(),
        selected = safe_grep_match(
          pattern = "alt",
          choices = gt_colnames()
        )
      )
      # Update marker ID select input allowing multiple selections
      updateSelectizeInput(
        server = TRUE,
        session = session,
        inputId = "marker_ID",
        choices = unique_marker_id()
      )
      # Update chromosome ID select input
      updateSelectizeInput(
        server = TRUE,
        session = session,
        inputId = "chr_ID",
        choices = unique_chrom()
      )
    })

    # DESIGN MARKERS --------------------------------------------------------

    # Designed markers and their alignment plots
    kasp_des.result <- reactiveVal(NULL)
    kasp_des.plot <- reactiveVal(NULL)

    # Design KASP Markers
    observeEvent(input$run_but, {
      data_source_field <- if (identical(input$upload_choice, "snpEff Annotated VCF")) {
        list("snpEff Annotated VCF File" = vcf_path())
      } else {
        list("Genotype Matrix File" = gt_data())
      }

      if (!validate_required_inputs(c(
        list(
          "Variant Data Type" = input$upload_choice,
          "Variant IDs Column" = input$variant_id_col,
          "Chromosome Column" = input$chrom_col,
          "Position Column" = input$pos_col,
          "Reference Allele Column" = input$ref_al_col,
          "Alternate Allele Column" = input$alt_al_col,
          "Genotype Data Start Column" = input$geno_start,
          "Marker Variants" = input$marker_ID,
          "Target Chromosome" = input$chr_ID,
          "Genome Reference File" = input$genome_file,
          "Minor Allele Frequency (MAF)" = input$maf,
          "Genome Version" = input$genome_version,
          "Trait" = input$trait,
          "Owner" = input$owner
        ),
        data_source_field
      ))) {
        return()
      }

      shinybusy::show_modal_spinner(
        spin = "fading-circle",
        color = "#27AE60",
        text = "Designing KASP Marker(s)... Please wait."
      )

      # Determine the data source based on user's choice
      vcf_file_path <- if (input$upload_choice == "snpEff Annotated VCF") {
        vcf_path()
      } else {
        NULL
      }

      genotype_df <- if (input$upload_choice == "Genotype Matrix (Processed)") {
        gt_data()
      } else {
        NULL
      }

      tryCatch({
        # Call the vectorized kasp_marker_design function once
        result_data <- kasp_marker_design(
          vcf_file = vcf_file_path,
          gt_df = genotype_df,
          variant_id_col = input$variant_id_col,
          chrom_col = input$chrom_col,
          pos_col = input$pos_col,
          ref_al_col = input$ref_al_col,
          alt_al_col = input$alt_al_col,
          geno_start = input$geno_start,
          marker_IDs = input$marker_ID, # Pass the vector of IDs
          chr = input$chr_ID,
          genome_file = input$genome_file$datapath,
          plot_file = tempdir(),
          save_alignment = FALSE,
          region_name = input$reg_name,
          maf = input$maf
        )

        # Generate Intertek table
        intertek_table <- make_intertek_table(
          marker_data = result_data,
          genome_version = input$genome_version,
          region_name = input$reg_name,
          trait = input$trait,
          owner = input$owner
        )

        # Store the consolidated results
        kasp_des.result(intertek_table)
        kasp_des.plot(result_data$plot)

        # Reveal the results
        shinyjs::hide("placeholder_ui")
        shinyjs::show("analysis_ui")
        bslib::nav_select("results_tabs", "table_tab")
        bslib::toggle_sidebar(id = "sidebar", open = "closed")

        shinyWidgets::show_alert(
          title = "Markers designed",
          text = paste(nrow(intertek_table), if (nrow(intertek_table) == 1) "KASP assay is" else "KASP assays are", "ready in the Marker Table."),
          type = "success",
          showCloseButton = TRUE,
          timer = 5000
        )
      }, error = function(e) {
        kasp_des.result(NULL)
        kasp_des.plot(NULL)
        shinyWidgets::show_alert(
          title = "Couldn't design the markers",
          text = paste0(e$message, " Check that the chromosome names in the variant file match the reference genome."),
          type = "error",
          showCloseButton = TRUE,
          timer = 8000
        )
      }, finally = {
        shinybusy::remove_modal_spinner()
      })
    })

    # MARKER TABLE ----------------------------------------------------------

    # Fill the alignment plot picker
    observe({
      req(kasp_des.result(), kasp_des.plot())

      updateSelectizeInput(session,
        inputId = "plot_choice",
        choices = names(kasp_des.plot())
      )
    })

    # Summary line above the table
    output$table_summary <- renderUI({
      req(kasp_des.result())
      n <- nrow(kasp_des.result())
      tags$p(
        class = "small text-muted",
        tags$b(n), if (n == 1) " assay" else " assays", " designed on ", tags$b(input$chr_ID),
        if (nzchar(input$reg_name)) tagList(" (region ", tags$b(input$reg_name), ")"),
        ". In the sequence, the target variant is shown in square brackets as [reference/alternate]."
      )
    })

    # Marker table
    output$kasp_table <- DT::renderDT({
      req(kasp_des.result())
      DT::datatable(kasp_des.result(),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })

    # Marker table download
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

    # ALIGNMENT PLOT --------------------------------------------------------

    # Plot for the selected marker
    output$plot <- renderPlot({
      req(kasp_des.plot(), input$plot_choice)
      print(kasp_des.plot()[[input$plot_choice]])
    })

    # All plots in one pdf
    output$download_plot <- downloadHandler(
      filename = function() {
        markers <- names(kasp_des.plot())
        tag <- if (length(markers) == 1) markers else paste0(length(markers), "_markers")
        paste0("alignment_", gsub("[^[:alnum:]_-]", "_", tag), ".pdf")
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
