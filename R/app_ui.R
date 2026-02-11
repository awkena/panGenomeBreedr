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

    # Floating card CSS
    tags$head(
      tags$style(HTML("
    /* The transition makes the movement smooth */
    .feature-card {
      transition: transform 0.3s ease-in-out, box-shadow 0.3s ease-in-out !important;
      cursor: pointer;
    }

    /* This creates the floating/lift effect */
    .feature-card:hover {
      transform: translateY(-10px);
      box-shadow: 0 10px 20px rgba(0,0,0,0.15) !important;
      border-color: #10b981 !important; /* Optional: highlight border on hover */
    }
  "))
    ),
    # Your application UI logic
    shinybusy::add_busy_bar(color = "#7f8c8d", height = "10px"),
    navbarPage(
      collapsible = TRUE,
      theme = bslib::bs_theme(
        version = 5,
        bootswatch = "flatly",
        primary = "#1e3a5f",
        success = "#27AE60",
        info = "#3498DB",
        warning = "#F39C12",
        danger = "#E74C3C"
        # primary = "#145214"  ,#"#228B22", # "#2E8B57",
        # success = "#FF8A65"
      ),
      fluid = TRUE,
      title = "panGenomeBreedr",
      id = "nav_bar",
        header = tagList(
          tags$head(
            # Helper script for redirecting users to vignette to learn more
            tags$script('
    Shiny.addCustomMessageHandler("open_link", function(message) {
    window.open(message.url, "_blank");
    });
    '),)
          ),

      # Home Tab
      nav_panel(
        title = tagList(icon("home"), "Home"),
        value = "home",

        # Welcome Section
        tags$section(
          class = "text-center p-5 bg-light border-bottom",
          tags$h2("Welcome to panGenomeBreedr!", class = "display-5 fw-bold mb-4"),
          tags$p(
            class = "lead text-muted mx-auto",
            style = "max-width: 800px;",
            "panGB is your unified platform for modern, crop-agnostic pangenome-enabled breeding.
            Engineered to simplify variant analysis, marker design, and validation."
          )
        ),

        # Conceptual Workflow Section (Embedded Image)
        tags$section(
          class = "text-center py-5 bg-white",
          tags$h3("Conceptual Workflow of panGB", class = "fw-bold mb-4"),
          tags$img(
            src = "www/workflow_ui.png", # Path relative to inst/app/www
            alt = "panGB Workflow Diagram",
            class = "img-fluid shadow-sm rounded-3",
            style = "max-width: 850px; border: 1px solid #eee;"
          )
        ),

        # 3. Feature Cards (Using bslib::layout_column_wrap)
        tags$section(
          class = "container py-5",
          bslib::layout_column_wrap(
            width = 1/4, # 4 columns
            gap = "20px",

            # Variant Discovery Card
            bslib::card(
              class = "feature-card text-center shadow-sm border-0 p-3",
              tags$div(class = "h1 text-primary mb-3", icon("magnifying-glass")),
              tags$h5("Variant Discovery", class = "fw-bold"),
              bslib::card_footer(
                actionButton("btn_variant_discovery", "Learn more", class = "btn-success rounded-pill px-4 btn-sm"),
                class = "bg-transparent border-0"
              )
            ),

            # KASP Marker Design Card
            bslib::card(
              class = "feature-card text-center shadow-sm border-0 p-3",
              tags$div(class = "h1 text-primary mb-3", icon("dna")),
              tags$h5("KASP Marker Design", class = "fw-bold"),
              bslib::card_footer(
                actionButton("btn_kasp_marker", "Learn more", class = "btn-success rounded-pill px-4 btn-sm"),
                class = "bg-transparent border-0"
              )
            ),

            # Marker Validation Card
            bslib::card(
              class = "feature-card text-center shadow-sm border-0 p-3",
              tags$div(class = "h1 text-primary mb-3", icon("check-circle")),
              tags$h5("KASP Marker Validation", class = "fw-bold"),
              bslib::card_footer(
                actionButton("btn_marker_validation", "Learn more", class = "btn-success rounded-pill px-4 btn-sm"),
                class = "bg-transparent border-0"
              )
            ),

            # Introgression Support card
            bslib::card(
              class = "feature-card text-center shadow-sm border-0 p-3",
              tags$div(class = "h1 text-primary mb-3", icon("cogs")),
              tags$h5("Introgression Support", class = "fw-bold"),
              bslib::card_footer(
                actionButton("btn_introgression", "Learn more", class = "btn-success rounded-pill px-4 btn-sm"),
                class = "bg-transparent border-0"
              )
            )
          )
        ),

        # Footer Section
        tags$footer(
          class = "py-5 border-top bg-light text-center",
          tags$p(class = "text-muted small mb-1", HTML("&copy; 2025 panGenomeBreedr")),
          tags$p(
            class = "text-muted extra-small mx-auto",
            style = "max-width: 700px; font-size: 0.8rem;",
            "Developed by: Alexander Wireko Kena, Israel Tawiah Tetteh, Cruet Burgos, Fanna Maina, Linly Banda, Jacques Faye, Benjamin Annor, Terry Felderhoff, Geoffrey Preston Norris"
          )
        )
      ),

      bslib::nav_item(), bslib::nav_item(), # space home tab from the rest of the tab

      ## Variant Annotation Tab
      bslib::nav_panel(
        title = "Variant Discovery",
        icon = icon("magnifying-glass"),
        mod_variant_discovery_ui("variant_discovery_1")
      ),
      bslib::nav_item(),

      ## KASP Marker Design UI
      bslib::nav_panel(
        title = "KASP Marker Design",
        icon = icon("dna"),
        mod_kasp_marker_design_ui("kasp_marker_design_1")
      ),
      bslib::nav_item(), # space KASP marker design tab from the rest of the tab

      ## Marker validation tab with sub functionalities.
      bslib::nav_panel(
        title = "KASP Marker Validation",
        icon = icon("check-circle"),
        bslib::navset_card_pill(
          id = "marker_validation_tabs",

          # Read Kasp file Tab
          bslib::nav_panel(
            title = tags$strong("Upload Data"),
            icon = icon("file-import"),
            bslib::card(
              height = "700px",
              bslib::card_header(h5(tags$b("KASP Genotyping File"))),
              mod_mv_read_kasp_csv_ui("mv_read_kasp_csv_1")
            )
          ),

          # Tab for Kasp file Summary / Statistics
          bslib::nav_panel(
            title = tags$strong("Data Summary"),
            icon = icon("table"),
            bslib::layout_column_wrap(
              width = 1,
              bslib::card(
                bslib::card_header(h5(tags$b("Plate Statistics"))),
                mod_mv_nsamples_plate_ui("mv_nsamples_plate_1") # plate count
              ),
              bslib::card(
                bslib::card_header(h5(tags$b("Allele Call Summary"))),
                mod_mv_get_alleles_ui("mv_get_alleles_1") # allele info
              )
            )
          ),

          # Tab for color coding
          bslib::nav_panel(
            title = tags$strong("Color Code Calls"),
            icon = icon("palette"),
            bslib::card(
              bslib::card_header(h5(tags$b("Color Code Genotypes"))),
              mod_mv_kasp_color_ui("mv_kasp_color_1")
            )
          ),

          # Tab for QC plots
          bslib::nav_panel(
            title = tags$strong("QC Plots & Layout"),
            icon = icon("chart-bar"),
            bslib::card(
              bslib::card_header(h5(tags$b("QC Plots, Predictions & Plate Layout"))),
              mod_mv_kasp_qc_ggplot_ui("mv_kasp_qc_ggplot_1")
            )
          ),

          # Tab for Predictions & Plots
          bslib::nav_panel(
            title = tags$strong("Plate Prediction QC"),
            icon = icon("chart-line"),
            bslib::card(
              bslib::card_header(h5(tags$b("QC Prediction Summary & Plots"))),
              mod_mv_pred_sum_stat_ui("mv_pred_sum_stat_1")
            )
          )
        )
      ), bslib::nav_item(), # space marker validation tab from the rest of the tab

      ## Trait Introgression Decision Suite
      navbarMenu(
        title = "Trait Introgression Decision Suite",
        icon = icon("cogs"),

        # Tab for Trait Introgression Hypothesis Testing
        bslib::nav_panel(
          title = "Trait Introgression Hypothesis Testing",
          mod_ds_trait_introg_hypothesis_test_ui("ds_trait_introg_hypothesis_test_1")
        ),

        # Tab for MABC Decision support
        bslib::nav_panel(
          title = "Decision Support for MABC",
          mod_ds_marker_ass_bac_ui("ds_marker_ass_bac_1")
        ),

        # Tab for Foreground Selection
        bslib::nav_panel(
          title = "Foreground Selection Support",
          mod_ds_foreground_select_ui("ds_foreground_select_1")
        )
      ),
      bslib::nav_spacer(), # Push power button to far rught

      bslib::nav_item(
        actionButton(
          inputId = "power_btn",
          label = NULL,
          icon = icon("power-off"),
          class = "btn-outline-danger border-0",
          style = "margin-top: 5px;" # Align with other nav items
        )
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
      app_title = "panGenomeBreedr"
    ),
    # Reload instead of grey scaale of death when app crashes.
    tags$script(HTML("
  $(document).on('shiny:disconnected', function(event) {
    // Instead of showing the grey screen, just reload the app
    location.reload();
  });
"))

    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()

  )

}
