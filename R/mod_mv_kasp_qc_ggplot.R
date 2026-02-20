#' QC Plot UI Function
#'
#' @description A shiny Module for generating KASP QC plots with predictions.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList uiOutput selectInput numericInput actionButton downloadButton plotOutput div textInput
#' @importFrom bslib card card_header card_body card_footer accordion accordion_panel
#'
mod_mv_kasp_qc_ggplot_ui <- function(id) {
  ns <- NS(id)

  tagList(
    bslib::card(
      bslib::card_header("QC Plot Configuration", class = "bg-dark text-white"),
      bslib::input_switch(id = ns("configure"), label = "Configure Plot", value = FALSE),
      conditionalPanel(
        condition = paste0('input["', ns("configure"), '"] == true'),
        tagList(
          bslib::card_body(
            fluidRow(
              # Card 1: Core Data Mapping
              column(
                width = 3,
                bslib::card(
                  class = "shadow p",
                  height = "100%",
                  bslib::card_header(tags$b("Core Data Mapping"),
                                     class = "bg-primary text-center", style = "font-size:18px;"),
                  bslib::card_body(
                    selectInput(
                      inputId = ns("well_id"),
                      label = "Grouping Column",
                      choices = NULL,
                      multiple = FALSE,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("fam_id"),
                      label = "FAM Coordinate Column",
                      choices = NULL,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("Hex_id"),
                      label = "HEX Coordinate Column",
                      choices = NULL,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("snp_id"),
                      label = "SNP ID Column",
                      choices = NULL,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("geno_call"),
                      label = "Genotype Call Column",
                      choices = NULL,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("group_id"),
                      label = "Positive Control Column",
                      choices = NULL,
                      width = "100%"
                    )
                  )
                )
              ),

              # Card 2: Quality Thresholds
              column(
                width = 3,
                bslib::card(
                  class = "shadow p",
                  height = "100%",
                  bslib::card_header(tags$b("Quality Thresholds"),
                                     class = "bg-danger text-center", style = "font-size:18px;"),
                  bslib::card_body(
                    textInput(
                      inputId = ns("group_unknown"),
                      label = "Unknown Genotype Label",
                      value = "?",
                      width = "100%"
                    ),
                    textInput(
                      inputId = ns("unused"),
                      label = "Unused Well Label",
                      value = "?",
                      width = "100%"
                    ),
                    textInput(
                      inputId = ns("blank"),
                      label = "NTC / Blank Label",
                      value = "NTC",
                      width = "100%"
                    ),
                    textInput(
                      inputId = ns("uncallable"),
                      label = "Uncallable Label",
                      value = "uncallable",
                      width = "100%"
                    ),
                    textInput(
                      inputId = ns("others"),
                      label = "Non-genotype Labels",
                      value = "Missing, Bad, Dupe, Over, Short",
                      width = "100%"
                    )
                  )
                )
              ),

              # Card 3: Visualization Parameters
              column(
                width = 3,
                bslib::card(
                  class = "shadow p",
                  height = "100%",
                  bslib::card_header(tags$b("Visualization Parameters"),
                                     class = "bg-info text-center",style = "font-size:18px;"),
                  bslib::card_body(
                    selectInput(
                      inputId = ns("scale"),
                      label = "Normalize Axes",
                      choices = c(TRUE, FALSE),
                      selected = TRUE,
                      width = "100%"
                    ),
                    numericInput(
                      inputId = ns("legendx_id"),
                      label = "Legend X Offset",
                      value = 0.6, min = 0, max = 1, step = 0.05,
                      width = "100%"
                    ),
                    numericInput(
                      inputId = ns("legendy_id"),
                      label = "Legend Y Offset",
                      value = 0.8, min = 0, max = 1, step = 0.05,
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("legend_box"),
                      label = "Legend Alignment",
                      choices = c("horizontal", "vertical"),
                      width = "100%"
                    ),
                    selectInput(
                      inputId = ns("legend_pos"),
                      label = "Legend Position",
                      choices = c("inside", "left", "right", "bottom", "top", "none"),
                      width = "100%"
                    )
                  )
                )
              ),

              # Card 4: Aesthetic Controls
              column(
                width = 3,
                bslib::card(
                  class = "shadow p",
                  height = "100%",
                  bslib::card_header(tags$b("Aesthetic Controls"),
                                     class = "bg-success text-center",
                                     style = "font-size:18px;"),
                  bslib::card_body(
                    div(
                      selectInput(
                        inputId = ns("pred_col_id"),
                        label = "Positive Control Colors",
                        choices = grDevices::colors(),
                        selected = c("black", "firebrick", "cornflowerblue", "beige"),
                        multiple = TRUE,
                        width = "100%"
                      ),
                      helpText(
                        tags$em("Required Order:"),
                        "Blank, False, True, Unverified"
                      )
                    ),
                    numericInput(
                      inputId = ns("textsize_id"),
                      label = "Plot Font Size",
                      value = 12, min = 1, max = 30, step = 1,
                      width = "100%"
                    ),
                    numericInput(
                      inputId = ns("alpha_id"),
                      label = "Color Opacity",
                      value = 0.8, min = 0, max = 1, step = 0.05,
                      width = "100%"
                    ),
                    numericInput(
                      inputId = ns("expand_id"),
                      label = "Axis Expansion",
                      value = 0.6, min = 0, max = 1, step = 0.05,
                      width = "100%"
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
    splitLayout(
      bslib::accordion(
        width = "100%",
        open = TRUE,
        bslib::accordion_panel(
          title = "KASP marker genotyping QC plot overlaid with predicitons",
          selectInput(
            inputId = ns("plate_choice"),
            label = "Plate Selection",
            choices = NULL,
            width = "60%"
          ),hr(),
          div(
            style = "display: flex; justify-content: center;",
            plotOutput(outputId = ns("qc_plot2"), width = "100%", height = "500px")
          ),
          card(card_footer(fluidRow(
            column(
              3,
              textInput(inputId = ns("file_name2"), label = "Enter Filename", value = "QC_plot_with_prediction")
            ),
            column(
              3,
              selectInput(
                inputId = ns("file_type2"), label = "File Format",
                choices = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff"),
                selected = "pdf"
              )
            ),
            column(
              3,
              numericInput(
                inputId = ns("width"),
                label = "Set Plot width",
                value = 6,
                min = 1
              )
            ),
            column(
              3,
              numericInput(
                inputId = ns("height"),
                label = "Set Plot height",
                value = 6,
                min = 1
              )
            )
          ), downloadButton(outputId = ns("download_plot2"), label = "Download Plot", class = "btn-success")))
        )
      ),
      bslib::accordion(
        open = TRUE,
        bslib::accordion_panel(
          title = "KASP Genotyping Plate Layout",
          br(),
          plotOutput(outputId = ns("plate_layout_plot"), width = "100%", height = "580px"),
          bslib::card(card_footer(
            fluidRow(
              column(
                3,
                textInput(
                  inputId = ns("file_name"),
                  label = "Enter Filename",
                  value = "Plate layout 1"
                )
              ),
              column(
                3,
                numericInput(
                  inputId = ns("width"),
                  label = "Set Plot Width",
                  value = 8, min = 1
                )
              ), column(
                3,
                numericInput(
                  inputId = ns("height"),
                  label = "Set Plot Height",
                  value = 5, min = 1
                )
              )
            ),
            downloadButton(
              outputId = ns("download_plateplot"),
              label = "Download Plot", class = "btn-success"
            )
          ))
        )
      )
    )
  )
}




#' mv_kasp_qc_ggplot Server Functions
#'
#' @param id An ID string that corresponds with the ID used to call the module's UI function
#' @param req_inputs A reactive that returns a list containing required input data (must contain drop_down, geno_call_id, and snp_id)
#' @param color_coded A reactive that returns the color-coded KASP data
#'
#' @return A reactive containing the generated QC plots
#'
#' @noRd
#'
#' @importFrom shiny moduleServer reactive reactiveVal observeEvent req updateSelectInput downloadHandler renderPlot
#' @importFrom ggplot2 ggsave
#' @importFrom stats setNames
#'
mod_mv_kasp_qc_ggplot_server <- function(id, kasp_data, color_coded) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # Populate drop down for FAM, HEX, SNPID & GENOTYPE CALLS
      names_pop <- reactive({
        req(kasp_data)
        colnames(kasp_data)
      })

      # Update select input.
      observe({
        req(names_pop(), color_coded())
        updateSelectInput(session,
          inputId = "group_id",
          choices = c("None", names_pop()),
          selected = "None"
        )
        updateSelectInput(session,
          inputId = "geno_call",
          choices = names_pop(),
          selected = safe_grep_match(
            pattern = "call",
            choices = names_pop()
          )
        )
        updateSelectInput(session,
          inputId = "snp_id",
          choices = names_pop(),
          selected = safe_grep_match(
            pattern = "snp",
            choices = names_pop()
          )
        )
        updateSelectInput(session,
          inputId = "well_id",
          choices = names_pop(),
          selected = safe_grep_match(
            pattern = "MasterWell",
            choices = names_pop()
          )
        )

        updateSelectInput(session,
          inputId = "Hex_id",
          choices = names_pop(),
          selected = safe_grep_match(
            pattern = "Y",
            choices = names_pop()
          )
        )
        updateSelectInput(session,
          inputId = "fam_id",
          choices = names_pop(),
          selected = safe_grep_match(
            pattern = "X",
            choices = names_pop()
          )
        )

        updateSelectInput(session,
          inputId = "plate_choice",
          choices = names(color_coded())
        )
      })


      result_plot <- reactiveVal(NULL)

      observe({
        # Validate all required inputs exist
        req(
          color_coded(),
          input$fam_id, input$Hex_id, input$geno_call,
          input$snp_id, input$group_id,
          input$legendx_id, input$legendy_id,
          input$legend_box, input$alpha_id,
          input$textsize_id, input$legend_pos,
          input$scale, input$pred_col_id,
          input$blank, input$uncallable,
          input$others, input$group_unknown,
          input$unused, length(input$pred_col_id) == 4
        )

        # Get the actual column names from the current data

        cc_res <- color_coded()
        req(length(cc_res) >= 1)

        # Get the actual column names from the current data
        current_cols <- colnames(cc_res[[1]])

        # Only proceed if selected inputs exist in the current data
        req(
          input$snp_id %in% current_cols,
          input$fam_id %in% current_cols,
          input$group_id %in% current_cols |input$group_id == "None",
          input$geno_call %in% current_cols,
          input$Hex_id %in% current_cols
        )

        tryCatch(
          {
            # Generate the plot
            result <- panGenomeBreedr::kasp_qc_ggplot2(
              pdf = FALSE,
              blank = input$blank,
              uncallable = input$uncallable,
              others = trimws(unlist(strsplit(input$others, ","))),
              Group_unknown = input$group_unknown,
              unused = input$unused,
              x = color_coded(),
              FAM = input$fam_id,
              HEX = input$Hex_id,
              geno_call = input$geno_call,
              snp_id = input$snp_id,
              Group_id = if (input$group_id == "None") NULL else input$group_id,
              legend.pos.x = input$legendx_id,
              legend.pos.y = input$legendy_id,
              legend.box = input$legend_box,
              alpha = input$alpha_id,
              text_size = input$textsize_id,
              legend.pos = input$legend_pos,
              scale = input$scale,
              pred_cols = c(
                 Blank = input$pred_col_id[1],
                 False = input$pred_col_id[2],
                 True = input$pred_col_id[3],
                 Unverified = input$pred_col_id[4]
              ),
              expand_axis = input$expand_id
            )
            result_plot(result)
          },
          error = function(e) {
            # Show error message
            shinyjs::delay(100, {
              shinyWidgets::show_alert(
                title = "Plot Generation Error",
                text = paste("Failed to generate plot:", e$message),
                type = "error",
                showCloseButton = TRUE,
                timer = 5000
              )
            })
            NULL
          }
        )
      })
      # Plate layout plot script.
      # Reactive for generating plate plot with error handling
      plot_plate_result <- reactive({
        # Validate inputs
        req(
          color_coded(), input$snp_id,
          input$geno_call, input$well_id,
          input$plate_choice
        )

        tryCatch(
          {
            # Generate plot
            plot_plate(
              x = color_coded()[input$plate_choice],
              well = input$well_id,
              geno_call = input$geno_call,
              snp_id = input$snp_id,
              filename = NULL,
              pdf = FALSE
            )
          },
          error = function(e) {
            shinyWidgets::show_alert(
              title = "Error!",
              text = paste("Failed to generate plot:", e$message),
              type = "error",
              showCloseButton = TRUE,
              timer = 5000
            )
            return(NULL)
          }
        )
      })

      # Render plot for plate layout.
      output$plate_layout_plot <- renderPlot({
        req(plot_plate_result(), input$plate_choice)

        plot_obj <- plot_plate_result()[[input$plate_choice]]

        if (is.null(plot_obj)) {
          shinyWidgets::show_toast(
            title = "",
            text = "No plot found for selected plate",
            type = "error"
          )
          return(NULL)
        }

        print(plot_obj)
      })

      # Render QC plot
      output$qc_plot2 <- renderPlot({
        req(result_plot(), input$plate_choice)
      tryCatch({
        print(result_plot()[input$plate_choice])

        shinyWidgets::show_toast(
          title = "",
          text = "Plot generated successfully",
          type = "success",
          timer = 2000,
          position = "bottom-end"
        )
      },error = function(e){
        shinyWidgets::show_toast(
          title = "",
          text = e$message,
          type = "error",
          timer = 2000,
          position = "bottom-end"
        )
      })
      })


      output$download_plot2 <- downloadHandler(
        filename = function() {
          clean_name <- gsub("[^[:alnum:]_-]", "_", input$file_name2)
          paste0(clean_name, ".pdf") # Force PDF for multi-page
        },
        content = function(file) {
          req(result_plot()) # Only require plot data (no plate_choice dependency)

          tryCatch(
            {
              # Start multi-page PDF
              grDevices::pdf(file,
                width = input$width,
                height = input$height,
                onefile = TRUE
              ) # Critical for multi-page

              # Print ALL plots in result_plot()
              for (plot in result_plot()) {
                print(plot)
              }
              # print(result_plot()[[input$plate_choice]])

              grDevices::dev.off()
            },
            error = function(e) {
              shinyWidgets::show_toast(
                paste("Download failed:", e$message),
                type = "error"
              )
            }
          )
        }
      )

      output$download_plateplot <- downloadHandler(
        filename = function() {
          req(input$file_name)
          paste0(input$file_name, ".pdf") # Single PDF containing ALL plots
        },
        content = function(file) {
          req(plot_plate_result()) # Ensure plots exist

          tryCatch(
            {
              # Start PDF (onefile=TRUE ensures multi-page)
              grDevices::pdf(file,
                width = input$width,
                height = input$height,
                onefile = TRUE
              )

              # Print ALL plots regardless of type
              for (plot_obj in plot_plate_result()) {
                # if (inherits(plot_obj, "grob")) {
                #   grid::grid.draw(plot_obj)  # Handle grid graphics
                # } else {
                print(plot_obj) # Handle ggplot2/base R plots
                # }
              }

              grDevices::dev.off()
            },
            error = function(e) {
              shinyWidgets::show_toast(
                title = "",
                type = "error",
                text = paste("Download failed:", e$message),
                timer = 5000
              )
            }
          )
        }
      )
    }
  )
}

## To be copied in the UI
# mod_mv_kasp_qc_ggplot_ui("mv_kasp_qc_ggplot_1")

## To be copied in the server
# mod_mv_kasp_qc_ggplot_server("mv_kasp_qc_ggplot_1")
