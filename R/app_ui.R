#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bslib bs_theme nav_panel nav_item card card_header card_footer navset_card_pill layout_column_wrap
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    navbarPage(
      theme = bslib::bs_theme(
        version = 5,
        bootswatch = "flatly",
        # heading_font = font_google("Open Sans", wght = 700),
        # code_font = font_google("Roboto", wght = 500),
        # base_font = font_google("Roboto"),
        # heading_font = font_google("Roboto Condensed"),
        primary = "#303F9F",
        success = "#FF8A65"
      ),
      fluid = TRUE,
      title = "PanGB prototype 1",
      id = "nav_bar",
      header = tags$div(
        style = "position: absolute; right: 10px; top: 20px;",
        actionButton(
          inputId = "refresh_btn",
          label = "Refresh",
          icon = icon("sync-alt"),
          style = "background-color: transparent; border: none;"
        )
      ),

      # Home Tab
      bslib::nav_panel(
        title = "Home",
        icon = icon("home")
      ),
      bslib::nav_item(), bslib::nav_item(), # space home tab from the rest of the tab
      # Variant Annotation Tab
      bslib::nav_panel(
        title = "Variant Discovery",
        icon = icon("microscope"),
        mod_variant_discovery_ui("variant_discovery_1")
      ),
      bslib::nav_item(),

      # KASP Marker Design UI
      bslib::nav_panel(
        title = "KASP Marker Design",
        icon = icon("dna"),
        mod_kasp_marker_design_ui(id = "kasp_marker_design_1") # kasp marker design UI module
      ),
      bslib::nav_item(),
      # Marker validation tab with sub functionalities.
      bslib::nav_panel(
        title = "Marker Validation",
        icon = icon("check-circle"),
        bslib::navset_card_pill(
          # Tab for kasp file import
          bslib::nav_panel(
            title = tags$strong("Import Data"),
            icon = icon("file-import"),
            bslib::card(
              bslib::card_header(h5(tags$b("Import KASP Genotyping File"))),
              mod_mv_read_kasp_csv_ui("mv_read_kasp_csv_1")
            )
          ),

          # Tab for data summary
          bslib::nav_panel(
            title = tags$strong("Data Summary"),
            icon = icon("table"),
            bslib::layout_column_wrap(
              width = 1,
              bslib::card(
                bslib::card_header(h5(tags$b("Sample Statistics"))),
                mod_mv_nsamples_plate_ui("mv_nsamples_plate_1")
              ),
              bslib::card(
                bslib::card_header(h5(tags$b("Allele Information"))),
                mod_mv_get_alleles_ui("mv_get_alleles_1")
              )
            )
          ),

          # Tab for color coding
          bslib::nav_panel(
            title = tags$strong("Color Code"),
            icon = icon("palette"),
            bslib::card(
              bslib::card_header(h5(tags$b("Color Code Calls"))),
              mod_mv_kasp_color_ui("mv_kasp_color_1")
            )
          ),
          # Tab for Visualization
          bslib::nav_panel(
            title = tags$strong("QC Plots"),
            icon = icon("chart-bar"),
            bslib::card(
              bslib::card_header(h5(tags$b("Quality Control Plots Overlayed with Predictions"))),
              mod_mv_kasp_qc_ggplot_ui("mv_kasp_qc_ggplot_1")
            )
          ),
          # Tab to generate predictions
          bslib::nav_panel(
            title = tags$strong("QC Predictions"),
            icon = icon("chart-line"),
            bslib::card(
              bslib::card_header(h5(tags$b("Quality Control Predictions"))),
              mod_mv_pred_sum_stat_ui("mv_pred_sum_stat_1")
            )
          ),
          bslib::nav_panel(
            title = tags$strong("Plate Layout"),
            icon = icon("chart-bar"),
            bslib::card(
              bslib::card_header(h5(tags$b("Plot Plate Design"))),
              mod_mv_plate_plot_ui("mv_plate_plot_1")
            )
          )
        )
      ), bslib::nav_item(),
      # # Other Breeder Centered functionalities
      navbarMenu(
        title = "Decision Support",
        icon = icon("cogs"),
        bslib::nav_panel(
          title = "Trait Introgression & MABC",
          navset_card_pill(
            nav_panel(
              title = "Import & Process Data",
              mod_ds_process_data_ui("ds_process_data_1")
            ),
            nav_panel(
              title = "Generate Heatmap",
              mod_ds_create_heatmap_ui(id = "ds_create_heatmap_1")
            ),
            nav_panel(title = "Introgression Test"),
            nav_panel(title = "Decision Support for MABC")
          )
        ),
        # mod for this could appear here
        bslib::nav_panel(
          title = "Foreground Section"
        )
      ),
      bslib::nav_item(),
      # More tab for other information
      navbarMenu(
        title = "More",
        icon = icon("ellipsis-v"),

        # Tab for issues report
        bslib::nav_panel(
          title = "Report Issues",
        ),
        # Help tab
        bslib::nav_panel(title = "Help")
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "PanGB.webapp"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
