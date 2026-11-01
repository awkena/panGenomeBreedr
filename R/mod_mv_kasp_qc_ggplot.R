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

  ui_sidebar_layout(

    # SIDEBAR ---------------------------------------------------------------
    sidebar = tagList(
      selectInput(
        inputId = ns("plate_choice"),
        label = "Plate",
        choices = NULL,
        width = "100%"
      ),
      bslib::accordion(
        open = "Columns",
        multiple = TRUE,

        # Column mapping
        bslib::accordion_panel(
          title = "Columns",
          icon = icon("table-columns"),
          selectInput(
            inputId = ns("well_id"),
            label = label_help("Well Column", "Column with the well position (e.g. A01), used for the plate layout."),
            choices = NULL,
            multiple = FALSE,
            width = "100%"
          ),
          selectInput(
            inputId = ns("fam_id"),
            label = label_help("FAM Coordinate Column", "Column with the FAM fluorescence values (x axis)."),
            choices = NULL,
            width = "100%"
          ),
          selectInput(
            inputId = ns("Hex_id"),
            label = label_help("HEX Coordinate Column", "Column with the HEX fluorescence values (y axis)."),
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
            label = label_help(
              "Positive Control Column",
              "Optional column with the expected genotype of positive controls. When set, the plot marks whether each control was predicted correctly."
            ),
            choices = NULL,
            width = "100%"
          )
        ),

        # Labels used in the call column
        bslib::accordion_panel(
          title = "Call labels",
          icon = icon("tags"),
          textInput(
            inputId = ns("group_unknown"),
            label = label_help("Unknown Genotype Label", "Label used in the positive control column when the expected genotype is unknown."),
            value = "?",
            width = "100%"
          ),
          textInput(
            inputId = ns("unused"),
            label = label_help("Unused Well Label", "Label used for empty or unused wells."),
            value = "?",
            width = "100%"
          ),
          textInput(
            inputId = ns("blank"),
            label = label_help("NTC / Blank Label", "Label of the no-template-control (blank) wells."),
            value = "NTC",
            width = "100%"
          ),
          textInput(
            inputId = ns("uncallable"),
            label = label_help("Uncallable Label", "Label used for wells the software could not call."),
            value = "Uncallable",
            width = "100%"
          ),
          textInput(
            inputId = ns("others"),
            label = label_help("Non-genotype Labels", "Other labels in the call column that are not genotypes, separated by commas."),
            value = "Missing, Bad, Dupe, Over, Short",
            width = "100%"
          )
        )
      )
    ),

    # MAIN PANEL ------------------------------------------------------------
    bslib::layout_columns(
      col_widths = bslib::breakpoints(sm = 12, xl = c(6, 6)),

      # QC plot card
      ui_result_card(
        "QC plot with predictions",
        full_screen = TRUE,
        actions = tagList(
          # Appearance popover
          plot_options_button(
            selectInput(
              inputId = ns("scale"),
              label = "Scale axes to 0-1",
              choices = c(TRUE, FALSE),
              selected = TRUE
            ),
            bslib::layout_columns(
              numericInput(ns("textsize_id"), "Font Size", value = 12, min = 1, max = 30, step = 1),
              numericInput(ns("alpha_id"), "Opacity", value = 0.8, min = 0, max = 1, step = 0.05)
            ),
            numericInput(
              inputId = ns("expand_id"),
              label = "Axis expansion (room for legend)",
              value = 0.6, min = 0, max = 1, step = 0.05
            ),
            bslib::layout_columns(
              selectInput(
                inputId = ns("legend_pos"),
                label = "Legend Position",
                choices = c("inside", "left", "right", "bottom", "top", "none")
              ),
              selectInput(
                inputId = ns("legend_box"),
                label = "Legend Layout",
                choices = c("horizontal", "vertical")
              )
            ),
            bslib::layout_columns(
              numericInput(ns("legendx_id"), "Legend X", value = 0.6, min = 0, max = 1, step = 0.05),
              numericInput(ns("legendy_id"), "Legend Y", value = 0.8, min = 0, max = 1, step = 0.05)
            ),
            tags$p(class = "fw-semibold small mb-1", "Positive control colours"),
            bslib::layout_columns(
              colourpicker::colourInput(ns("col_blank"), "Blank", value = "black", returnName = TRUE),
              colourpicker::colourInput(ns("col_false"), "False", value = "firebrick", returnName = TRUE)
            ),
            bslib::layout_columns(
              colourpicker::colourInput(ns("col_true"), "True", value = "cornflowerblue", returnName = TRUE),
              colourpicker::colourInput(ns("col_unverified"), "Unverified", value = "beige", returnName = TRUE)
            )
          ),
          plot_download_button(
            ns, "download_plot2",
            file_name = "QC_plot_with_prediction",
            formats = c("PDF (all plates)" = "pdf", "PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff"),
            width = 6, height = 6,
            note = "PDF saves every plate; image formats save the selected plate."
          )
        ),
        # Plot with hover tooltip on top
        div(
          style = "position: relative; display: flex; justify-content: center;",
          plotOutput(
            outputId = ns("qc_plot2"),
            width = "100%",
            height = "500px",
            hover = hoverOpts(id = ns("plot_hover"), delay = 100, delayType = "debounce")
          ),
          uiOutput(ns("hover_tooltip"))
        )
      ),

      # Plate layout card
      ui_result_card(
        "Plate layout",
        full_screen = TRUE,
        actions = plot_download_button(
          ns, "download_plateplot",
          file_name = "Plate_layout",
          formats = c("PDF (all plates)" = "pdf"),
          width = 8, height = 5
        ),
        plotOutput(outputId = ns("plate_layout_plot"), width = "100%", height = "500px")
      )
    )
  )
}


#' mv_kasp_qc_ggplot Server Functions
#'
#' @param id An ID string that corresponds with the ID used to call the module's UI function
#' @param kasp_data A reactive that returns the imported KASP data
#' @param color_coded A reactive that returns the color-coded KASP data
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

      # COLUMN DROPDOWNS ------------------------------------------------------

      # Column names of the KASP data
      names_pop <- reactive({
        req(kasp_data())
        colnames(kasp_data())
      })

      # Fill the dropdowns and guess the right column for each
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
          selected = safe_grep_match(pattern = "call", choices = names_pop())
        )
        updateSelectInput(session,
          inputId = "snp_id",
          choices = names_pop(),
          selected = safe_grep_match(pattern = "snp", choices = names_pop())
        )
        updateSelectInput(session,
          inputId = "well_id",
          choices = names_pop(),
          selected = safe_grep_match(pattern = "MasterWell", choices = names_pop())
        )
        updateSelectInput(session,
          inputId = "Hex_id",
          choices = names_pop(),
          selected = safe_grep_match(pattern = "Y", choices = names_pop())
        )
        updateSelectInput(session,
          inputId = "fam_id",
          choices = names_pop(),
          selected = safe_grep_match(pattern = "X", choices = names_pop())
        )
        updateSelectInput(session,
          inputId = "plate_choice",
          choices = names(color_coded())
        )
      })

      # QC PLOT ---------------------------------------------------------------

      result_plot <- reactiveVal(NULL)

      # Build QC plots for all plates
      observe({
        req(
          color_coded(),
          input$fam_id, input$Hex_id, input$geno_call,
          input$snp_id, input$group_id,
          input$legendx_id, input$legendy_id,
          input$legend_box, input$alpha_id,
          input$textsize_id, input$legend_pos,
          input$scale,
          input$col_blank, input$col_false, input$col_true, input$col_unverified,
          input$blank, input$uncallable,
          input$others, input$group_unknown,
          input$unused
        )

        cc_res <- color_coded()
        req(length(cc_res) >= 1)

        # Only go on if the chosen columns are in the data
        current_cols <- colnames(cc_res[[1]])
        req(
          input$snp_id %in% current_cols,
          input$fam_id %in% current_cols,
          input$group_id %in% current_cols | input$group_id == "None",
          input$geno_call %in% current_cols,
          input$Hex_id %in% current_cols
        )

        tryCatch(
          {
            result <- kasp_qc_ggplot2(
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
              scale = as.logical(input$scale),
              pred_cols = c(
                Blank = input$col_blank,
                False = input$col_false,
                True = input$col_true,
                Unverified = input$col_unverified
              ),
              expand_axis = input$expand_id
            )
            result_plot(result)
          },
          error = function(e) {
            shinyjs::delay(100, {
              shinyWidgets::show_alert(
                title = "Couldn't draw the QC plot",
                text = e$message,
                type = "error",
                showCloseButton = TRUE,
                timer = 5000
              )
            })
            NULL
          }
        )
      })

      # Show plot for the selected plate
      output$qc_plot2 <- renderPlot({
        validate(need(kasp_data(), msg_no_data))
        validate(need(if_ready(color_coded()), msg_no_color))
        req(result_plot(), input$plate_choice)
        tryCatch({
          result_plot()[[input$plate_choice]]
        }, error = function(e) {
          shinyWidgets::show_toast(
            title = "",
            text = paste("Couldn't show the QC plot:", e$message),
            type = "error",
            timer = 4000,
            position = "bottom-end"
          )
          NULL
        })
      })

      # Hover tooltip with the call and group of the nearest point
      output$hover_tooltip <- renderUI({
        hover <- input$plot_hover
        req(hover, result_plot(), input$plate_choice)

        plt <- result_plot()[[input$plate_choice]]
        plot_data <- plt$data

        point <- nearPoints(plot_data, hover, xvar = input$fam_id, yvar = input$Hex_id, threshold = 10, maxpoints = 1)
        if (nrow(point) == 0) return(NULL)

        # Place it just off the mouse
        left_px <- hover$coords_css$x + 15
        top_px <- hover$coords_css$y + 15

        style <- paste0(
          "position: absolute; z-index: 1000; background-color: rgba(255, 255, 255, 0.95); ",
          "padding: 8px 12px; border: 1px solid #ccc; border-radius: 5px; ",
          "box-shadow: 2px 2px 10px rgba(0,0,0,0.2); pointer-events: none; ",
          "font-size: 13px; left: ", left_px, "px; top: ", top_px, "px;"
        )

        # Use tooltip_text from the plot data if it's there
        tooltip_str <- if ("tooltip_text" %in% colnames(point)) {
          point$tooltip_text
        } else {
          paste0("Sample: ", point[[input$well_id]], "\nCall: ", point[[input$geno_call]])
        }

        div(
          style = style,
          tags$span(style = "white-space: pre-wrap; font-weight: bold;", tooltip_str)
        )
      })

      # PLATE LAYOUT ----------------------------------------------------------

      # Build plate layouts. Returns the error message as text if it fails
      plot_plate_result <- reactive({
        req(
          color_coded(),
          input$snp_id,
          input$geno_call,
          input$well_id,
          input$plate_choice
        )

        tryCatch(
          {
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
            return(e$message)
          }
        )
      })

      # Show layout for the selected plate, or a ghost plate with the error
      output$plate_layout_plot <- shiny::renderPlot({
        validate(need(kasp_data(), msg_no_data))
        validate(need(if_ready(color_coded()), msg_no_color))
        res <- plot_plate_result()

        if (is.character(res)) {
          p_err <- render_ghost_plate_error(error_msg = res, text_size = 12)
          return(print(p_err))
        }

        req(input$plate_choice)

        plot_name <- grep(
          pattern = input$plate_choice,
          x = names(res),
          ignore.case = TRUE,
          value = TRUE
        )

        if (length(plot_name) == 0) {
          p_miss <- render_ghost_plate_error(
            "No matching plate name found in data",
            text_size = 12
          )
          return(print(p_miss))
        }

        plot_obj <- res[[plot_name]]

        if (is.null(plot_obj)) {
          p_null <- render_ghost_plate_error("Plot object was NULL", text_size = 12)
          return(print(p_null))
        }

        print(plot_obj)
      })

      # DOWNLOADS -------------------------------------------------------------

      # QC plot: PDF has all plates, images have the selected plate
      output$download_plot2 <- downloadHandler(
        filename = function() {
          clean_name <- gsub("[^[:alnum:]_-]", "_", input$download_plot2_name)
          ext <- input$download_plot2_type
          paste0(clean_name, ".", ext)
        },
        content = function(file) {
          req(result_plot())

          tryCatch(
            {
              ext <- input$download_plot2_type
              if (ext == "pdf") {
                grDevices::pdf(file,
                  width = input$download_plot2_width,
                  height = input$download_plot2_height,
                  onefile = TRUE
                )
                for (plot in result_plot()) { print(plot) }
                grDevices::dev.off()
              } else {
                req(input$plate_choice)
                current_plot <- result_plot()[[input$plate_choice]]
                ggplot2::ggsave(filename = file, plot = current_plot, device = ext,
                                width = input$download_plot2_width,
                                height = input$download_plot2_height)
              }
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

      # Plate layout: one PDF with every plate
      output$download_plateplot <- downloadHandler(
        filename = function() {
          req(input$download_plateplot_name)
          paste0(input$download_plateplot_name, ".pdf")
        },
        content = function(file) {
          req(plot_plate_result())

          tryCatch(
            {
              grDevices::pdf(file,
                width = input$download_plateplot_width,
                height = input$download_plateplot_height,
                onefile = TRUE
              )
              for (plot_obj in plot_plate_result()) {
                print(plot_obj)
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
