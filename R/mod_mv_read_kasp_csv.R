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
#' @importFrom shiny NS tagList sidebarLayout sidebarPanel mainPanel fileInput radioButtons numericInput textInput div actionButton icon
#' @importFrom bslib accordion accordion_panel navset_card_tab nav_panel
#'
mod_mv_read_kasp_csv_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        # Widget to upload kasp csv file
        fileInput(ns("Kasp_csv.file"),
          label = "Upload KASP Genotyping Results",
          accept = ".csv"
        ),
        # Widget to choose data type
        radioButtons(
          inputId = ns("datatype"),
          label = "Choose Data Format",
          choices = c("raw", "polished"),
          selected = "raw"
        ),
        # Widget for spacing
        numericInput(ns("data_space"),
          label = "Rows Between Segments",
          value = 2, min = 0, max = 100, step = 1
        ),
        # Widget for row tags
        textInput(
          inputId = ns("row_tags"),
          label = "Row Tags (ordered, comma-separated)",
          value = "Statistics , DNA , SNPs , Scaling , Data"
        ),
        # action button widget
        div(
          style = "display: flex; justify-content: center;",
          actionButton(ns("submit_btn"),
            label = "Submit",
            icon = icon("rocket"), class = "btn-info", width = "80%"
          )
        )
      ),
      mainPanel(
        # Display results from after reading csv.
        bslib::accordion(height = '500px',
          open = TRUE, # Default open section
          bslib::accordion_panel(title =   "Uploaded File Preview",
            navset_card_tab(height = '500px',
              nav_panel(title = "Data", DT::DTOutput(ns("kasp_data"),height = '500px')),
              nav_panel(title = "Statistics", DT::DTOutput(ns("kasp_statistics"))),
              nav_panel(title = "SNPS", DT::DTOutput(ns("kasp_snps"))),
              nav_panel(title = "DNA", DT::DTOutput(ns("kasp_DNA"))),
              nav_panel(title = "Scaling", DT::DTOutput(ns("kasp_scaling")))
            )
          )
        )
      )
    )
  )
}



#' mv_read_kasp_csv Server Functions
#' @import shiny
#'
#' @noRd

mod_mv_read_kasp_csv_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Create a null reactive component to store results from the read_kasp_csv()
    import_data <- reactiveVal(NULL)


    # get the list data from the uploaded KASP genotyping file.
    observeEvent(input$submit_btn, {
      req(input$row_tags, input$Kasp_csv.file, input$data_space, input$datatype)
      tryCatch(
        {
          # Split row tags and trim whitespace
          row_tags_list <- trimws(unlist(strsplit(x = input$row_tags, split = ",")))

          # Process the CSV
          read_kasp_result <- read_kasp_csv(
            file = input$Kasp_csv.file$datapath,
            spacing = input$data_space,
            data_type = input$datatype,
            row_tags = row_tags_list
          )

          # Update import data
          import_data(read_kasp_result)
          shinyWidgets::show_toast(
            title = "Success",
            text = "Data Imported Successfully",
            type = "success",
            timer = 2000,
            position = "bottom-end"
          )
        },
        error = function(e) {
          shinyWidgets::show_alert(
            title = "Error",
            text = paste("Error processing KASP data:", e$message),
            type = "error",
            showCloseButton = TRUE
          )
        }
      )
    })


    # Add plates if absent.
    add_plates <- reactiveVal(NULL)

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
            title = "Error",
            text = paste("Error processing plate data:", e$message),
            type = "error",
            showCloseButton = TRUE
          )
        }
      )
    })

    # Data ouptput
    output$kasp_data <- DT::renderDT({
      req(add_plates())
      tryCatch(
        DT::datatable(add_plates(),
          options = list(pageLength = 10, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "Data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })


    # Statistics output
    output$kasp_statistics <- DT::renderDT({
      req(import_data())
      tryCatch(
        DT::datatable(import_data()$Statistics,
          options = list(pageLength = 10, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "Statistics not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # SNPS output
    output$kasp_snps <- DT::renderDT({
      req(import_data())
      tryCatch(
        DT::datatable(import_data()$SNPs,
          options = list(pageLength = 10, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "SNPs data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # DNA output.
    output$kasp_DNA <- DT::renderDT({
      req(import_data())
      tryCatch(
        DT::datatable(import_data()$DNA,
          options = list(pageLength = 10, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "DNA data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # Scaling output
    output$kasp_scaling <- DT::renderDT({
      req(import_data())
      tryCatch(
        DT::datatable(import_data()$Scaling,
          options = list(pageLength = 10, scrollX = TRUE)
        ),
        error = function(e) {
          return(DT::datatable(data.frame(Message = "Scaling data not available"),
            options = list(dom = "t")
          ))
        }
      )
    })

    # Return reactive containing the data frame with plates added
    return(add_plates)
  })
}

## To be copied in the UI
# mod_mv_read_kasp_csv_ui("mv_read_kasp_csv_1")

## To be copied in the server
# mod_mv_read_kasp_csv_server("mv_read_kasp_csv_1")
