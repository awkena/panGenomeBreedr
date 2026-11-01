#' KASP CSV Data Reader UI Function
#'
#' @description A Shiny module for uploading and previewing KASP genotyping data in CSV format
#'
#' @param id Internal ID parameter for {shiny}
#'
#' @return A reactive object of Data to be used by other modules
#'
#' @noRd
#'
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel fileInput radioButtons numericInput textInput div actionButton icon uiOutput
#' @importFrom bslib accordion accordion_panel navset_card_tab nav_panel
#'
mod_mv_read_kasp_csv_ui <- function(id) {
  ns <- NS(id)
  ui_sidebar_layout(

    # SIDEBAR ---------------------------------------------------------------

    sidebar = tagList(
      # Upload KASP csv file, or use the example
      fileInput(ns("Kasp_csv.file"),
        label = "KASP Results File",
        accept = ".csv"
      ),
      example_data_link(ns("load_example")),
      uiOutput(ns("example_badge")),

      # Raw or polished file
      radioButtons(
        inputId = ns("datatype"),
        label = label_help(
          "Select Data Format",
          "raw: the file as exported by the genotyping service, with Statistics, DNA, SNPs, Scaling and Data sections. polished: a plain table with one row per well."
        ),
        choices = c("raw", "polished"),
        selected = "raw"
      ),

      # Raw files only: how the sections are laid out
      conditionalPanel(
        condition = sprintf("input['%s'] == 'raw'", ns("datatype")),
        bslib::accordion(
          open = FALSE,
          bslib::accordion_panel(
            title = "Raw file layout",
            icon = icon("table-list"),
            textInput(
              inputId = ns("row_tags"),
              label = label_help(
                "Ordered Row Tags",
                "Names of the sections in the raw file, in the order they appear, separated by commas."
              ),
              value = "Statistics , DNA , SNPs , Scaling , Data"
            ),
            numericInput(ns("data_space"),
              label = label_help(
                "Row Spacing",
                "Number of empty rows between two sections in the raw file."
              ),
              value = 2, min = 0, max = 100, step = 1
            )
          )
        )
      ),
      run_button(ns("submit_btn"), "Run Analysis")
    ),

    # MAIN PANEL ------------------------------------------------------------

    # One tab per section of the file
    bslib::navset_card_tab(
      title = tags$span(class = "fw-semibold", "Uploaded File Preview"),
      full_screen = TRUE,
      nav_panel(title = "Data", DT::DTOutput(ns("kasp_data"))),
      nav_panel(title = "Statistics", DT::DTOutput(ns("kasp_statistics"))),
      nav_panel(title = "SNPs", DT::DTOutput(ns("kasp_snps"))),
      nav_panel(title = "DNA", DT::DTOutput(ns("kasp_DNA"))),
      nav_panel(title = "Scaling", DT::DTOutput(ns("kasp_scaling")))
    )
  )
}


#' mv_read_kasp_csv Server Functions
#' @param demo A reactive that loads the example data each time it changes
#'   (used by the Home page "Try the demo" button).
#' @import shiny
#'
#' @noRd

mod_mv_read_kasp_csv_server <- function(id, demo = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # STATE -----------------------------------------------------------------

    # Result of read_kasp_csv()
    import_data <- reactiveVal(NULL)

    # Data with the plates column added. This is what the module returns
    add_plates <- reactiveVal(NULL)

    # File to read: the latest upload, or the bundled example
    data_path <- reactiveVal(NULL)
    using_example <- reactiveVal(FALSE)
    example_name <- "Genotyping_141.010_01.csv"
    default_row_tags <- "Statistics , DNA , SNPs , Scaling , Data"

    # READ FILE -------------------------------------------------------------

    # New upload
    observeEvent(input$Kasp_csv.file, {
      data_path(input$Kasp_csv.file$datapath)
      using_example(FALSE)
    })

    # Read a KASP file and store the result
    read_file <- function(path, datatype, row_tags, spacing) {
      tryCatch(
        {
          read_kasp_result <- if (datatype == "raw") {
            read_kasp_csv(
              file = path,
              row_tags = trimws(unlist(strsplit(row_tags, split = ","))),
              spacing = spacing,
              data_type = "raw"
            )
          } else {
            read_kasp_csv(file = path, data_type = "polished")
          }

          import_data(read_kasp_result)
          shinyWidgets::show_toast(
            title = "",
            text = "KASP data loaded",
            type = "success",
            timer = 2000,
            position = "bottom-end"
          )
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Couldn't read the KASP file",
            text = paste0(e$message, " Check that it is a KASP results file and that the data type matches."),
            type = "error",
            showCloseButton = TRUE
          )
        }
      )
    }

    # Run Analysis
    observeEvent(input$submit_btn, {
      if (!validate_required_inputs(list(
        "KASP Results File" = data_path(),
        "Data Format" = input$datatype,
        "Ordered Row Tags" = input$row_tags,
        "Row Spacing" = input$data_space
      ))) {
        return()
      }
      read_file(data_path(), input$datatype, input$row_tags, input$data_space)
    })

    # EXAMPLE DATA ----------------------------------------------------------

    # Load the vignette KASP file and run it straight away.
    # Called by the link on this tab or by `demo` (the Home page demo button)
    load_example <- function() {
      updateRadioButtons(session, "datatype", selected = "raw")
      updateTextInput(session, "row_tags", value = default_row_tags)
      updateNumericInput(session, "data_space", value = 2)

      data_path(example_file(example_name))
      using_example(TRUE)
      read_file(data_path(), "raw", default_row_tags, 2)
    }
    observeEvent(input$load_example, load_example())
    observeEvent(demo(), load_example(), ignoreInit = TRUE)

    # Clear the example
    observeEvent(input$clear_example, {
      data_path(NULL)
      using_example(FALSE)
      import_data(NULL)
      add_plates(NULL)
    })

    # "Using example data" badge
    output$example_badge <- renderUI({
      req(using_example())
      example_data_badge(example_name, clear_id = ns("clear_example"))
    })

    # ADD PLATES COLUMN -----------------------------------------------------

    # Add plates if absent
    observe({
      req(import_data())
      tryCatch(
        {
          if (is.list(import_data())) {
            refined_wplates <- plates_col(import_data()$Data) |> as.data.frame()
            add_plates(refined_wplates)
          } else {
            refined_wplates <- plates_col(import_data())
            add_plates(refined_wplates)
          }
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Couldn't read the plate layout",
            text = e$message,
            type = "error",
            showCloseButton = TRUE
          )
        }
      )
    })

    # PREVIEW TABLES --------------------------------------------------------

    # Data
    output$kasp_data <- DT::renderDT({
      validate(need(add_plates(), "Choose a KASP results file and click Run Analysis, or load the example data."))
      tryCatch(
        DT::datatable(add_plates(),
          options = list(pageLength = 5, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "Data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # Statistics
    output$kasp_statistics <- DT::renderDT({
      validate(need(import_data(), "Nothing to show yet."))
      tryCatch(
        DT::datatable(import_data()$Statistics,
          options = list(pageLength = 5, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "Statistics not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # SNPs
    output$kasp_snps <- DT::renderDT({
      validate(need(import_data(), "Nothing to show yet."))
      tryCatch(
        DT::datatable(import_data()$SNPs,
          options = list(pageLength = 5, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "SNPs data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # DNA
    output$kasp_DNA <- DT::renderDT({
      validate(need(import_data(), "Nothing to show yet."))
      tryCatch(
        DT::datatable(import_data()$DNA,
          options = list(pageLength = 5, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "DNA data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # Scaling
    output$kasp_scaling <- DT::renderDT({
      validate(need(import_data(), "Nothing to show yet."))
      tryCatch(
        DT::datatable(import_data()$Scaling,
          options = list(pageLength = 5, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "Scaling data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # Data with plates goes to the other Marker Validation modules
    return(add_plates)
  })
}

## To be copied in the UI
# mod_mv_read_kasp_csv_ui("mv_read_kasp_csv_1")

## To be copied in the server
# mod_mv_read_kasp_csv_server("mv_read_kasp_csv_1")
