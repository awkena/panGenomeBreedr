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
      header =
        tagList(
          tags$div(
            style = "position: absolute; right: 10px; top: 20px;",
            actionButton(
              inputId = "refresh_btn",
              label = "Refresh",
              icon = icon("sync-alt"),
              style = "background-color: transparent; border: none;"
            )
          ),
          tags$head(
            # Helper script for redirecting users to vignette to learn more
            tags$script('
  Shiny.addCustomMessageHandler("open_link", function(message) {
    window.open(message.url, "_blank");
  });
'),

            # Home page interface
            tags$style(HTML("
      .welcome-section {
        text-align: center;
        padding: 3rem 2rem;
        background: linear-gradient(to bottom, #f8f9fa, #ffffff);
      }
      .logo-hexagon {
        width: 80px;
        height: 80px;
        margin: 0 auto 1.5rem;
      }
      .welcome-title {
        font-size: 2rem;
        font-weight: 700;
        color: #1a202c;
        margin-bottom: 1.5rem;
      }
      .welcome-text {
        font-size: 1.1rem;
        color: #4a5568;
        max-width: 900px;
        margin: 0 auto;
        line-height: 1.6;
      }
      .workflow-section {
        padding: 3rem 2rem;
        background-color: #ffffff;
      }
      .workflow-title {
        text-align: center;
        font-size: 1.8rem;
        font-weight: 700;
        margin-bottom: 2rem;
        color: #1a202c;
      }
      .workflow-diagram {
        max-width: 100%;
        margin: 0 auto;
        display: block;
      }
      .card-container {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
        gap: 2rem;
        padding: 3rem 2rem;
        max-width: 1200px;
        margin: 0 auto;
      }
      .feature-card {
        background: #f7fafc;
        border-radius: 12px;
        padding: 2rem;
        text-align: center;
        transition: transform 0.3s ease, box-shadow 0.3s ease;
        border: 1px solid #e2e8f0;
      }
      .feature-card:hover {
        transform: translateY(-8px);
        box-shadow: 0 12px 24px rgba(0,0,0,0.15);
      }
      .feature-icon {
        width: 60px;
        height: 60px;
        margin: 0 auto 1rem;
      }
      .feature-title {
        font-size: 1.1rem;
        font-weight: 600;
        color: #2d3748;
        margin-bottom: 1rem;
      }
      .learn-more-btn {
        background-color: #10b981;
        color: white;
        border: none;
        border-radius: 20px;
        padding: 0.5rem 2rem;
        font-weight: 500;
        cursor: pointer;
        transition: background-color 0.3s ease;
      }
      .learn-more-btn:hover {
        background-color: #059669;
      }
      .why-section {
        padding: 3rem 2rem;
        background-color: #f8f9fa;
      }
      .why-title {
        text-align: center;
        font-size: 1.8rem;
        font-weight: 700;
        margin-bottom: 2rem;
        color: #1a202c;
      }
      .benefits-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
        gap: 2rem;
        max-width: 1000px;
        margin: 0 auto;
      }
      .benefit-card {
        background: white;
        padding: 1.5rem;
        border-left: 4px solid #10b981;
        border-radius: 4px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .benefit-title {
        font-weight: 700;
        color: #1a202c;
        margin-bottom: 0.5rem;
      }
      .benefit-text {
        color: #4a5568;
      }
      .footer-section {
        background-color: #ffffff;
        padding: 3rem 2rem;
        border-top: 1px solid #e2e8f0;
      }
      .footer-logos {
        text-align: center;
      }
      .footer-heading {
        font-size: 0.9rem;
        color: #718096;
        margin-bottom: 0.5rem;
        text-align: center;
      }
      .footer-org {
        font-size: 1.5rem;
        font-weight: 700;
        margin-bottom: 2rem;
        text-align: center;
      }
      .logo-row {
        display: flex;
        justify-content: center;
        align-items: center;
        flex-wrap: wrap;
        gap: 2rem;
        margin-bottom: 2rem;
      }
      .partner-logo {
        height: 60px;
        object-fit: contain;
      }
      .copyright {
        text-align: center;
        color: #718096;
        font-size: 0.9rem;
        margin-top: 2rem;
      }
      .developers {
        text-align: center;
        color: #718096;
        font-size: 0.85rem;
        margin-top: 0.5rem;
      }
    "))
          )
        ),

      # Home Tab
      nav_panel(
        title = tagList(icon("home"), "Home"),
        value = "home",

        # Welcome Section
        div(
          class = "welcome-section",
          # Logo placeholder - add image path in www folder
          # img(src = "logo_hexagon.png", class = "logo-hexagon"),
          h1(class = "welcome-title", "Welcome to panGenomeBreedr !"),
          p(
            class = "welcome-text",
            "panGB is your unified platform for modern, crop-agnostic pangenome-enabled breeding.
        Engineered to simplify natural and causal variant analysis, marker design, and marker
        validation, panGB empowers plant breeders and geneticists to unlock the full potential
        of pangenome resources in cultivar development."
          )
        ),

        # Workflow Section
        div(
          class = "workflow-section",
          h2(class = "workflow-title", "Conceptual Workflow of panGB"),
          # Workflow diagram
          # img(src = "www/workflow_diagram.png", class = "workflow-diagram")
        ),

        # Feature Cards
        div(
          class = "card-container",

          # Variant Discovery Card
          div(
            class = "feature-card",
            # Icon placeholder
            # img(src = "icon_variant_discovery.png", class = "feature-icon"),
            div(
              class = "feature-icon",
              style = "background-color: #e2e8f0; border-radius: 50%;
                   display: flex; align-items: center; justify-content: center;",
              icon("dna", style = "font-size: 30px; color: #718096;")
            ),
            div(class = "feature-title", "Variant Discovery"),
            actionButton(
              inputId = "btn_variant_discovery",
              label = "Learn more",
              class = "learn-more-btn"
            )
          ),

          # KASP Marker Design Card
          div(
            class = "feature-card",
            # Icon placeholder
            # img(src = "icon_kasp_marker.png", class = "feature-icon"),
            div(
              class = "feature-icon",
              style = "background-color: #e2e8f0; border-radius: 50%;
                   display: flex; align-items: center; justify-content: center;",
              icon("flask", style = "font-size: 30px; color: #718096;")
            ),
            div(class = "feature-title", "KASP Marker Design"),
            actionButton(
              "btn_kasp_marker",
              "Learn more",
              class = "learn-more-btn"
            )
          ),

          # Marker Validation Card
          div(
            class = "feature-card",
            # Icon placeholder
            # img(src = "icon_marker_validation.png", class = "feature-icon"),
            div(
              class = "feature-icon",
              style = "background-color: #e2e8f0; border-radius: 50%;
                   display: flex; align-items: center; justify-content: center;",
              icon("check-circle", style = "font-size: 30px; color: #718096;")
            ),
            div(class = "feature-title", "Marker Validation"),
            actionButton(
              "btn_marker_validation",
              "Learn more",
              class = "learn-more-btn"
            )
          ),

          # Introgression Support Card
          div(
            class = "feature-card",
            # Icon placeholder
            # img(src = "icon_introgression.png", class = "feature-icon"),
            div(
              class = "feature-icon",
              style = "background-color: #e2e8f0; border-radius: 50%;
                   display: flex; align-items: center; justify-content: center;",
              icon("cogs", style = "font-size: 30px; color: #718096;")
            ),
            div(class = "feature-title", "Introgression Support"),
            actionButton(
              "btn_introgression",
              "Learn more",
              class = "learn-more-btn"
            )
          )
        ),

        # Why panGB Section
        div(
          class = "why-section",
          h2(class = "why-title", "Why panGB?"),
          div(
            class = "benefits-grid",
            div(
              class = "benefit-card",
              div(class = "benefit-title", HTML("&#10003; Crop-agnostic platform")),
              div(class = "benefit-text", "for diverse breeding programs")
            ),
            div(
              class = "benefit-card",
              div(class = "benefit-title", HTML("&#10003; Platform-independent")),
              div(class = "benefit-text", "(Windows, MacOS)")
            ),
            div(
              class = "benefit-card",
              div(class = "benefit-title", HTML("&#10003; Free and Open Source")),
              div(class = "benefit-text", "Licensed under GNU GPL v3")
            )
            # div(
            #   class = "benefit-card",
            #   div(class = "benefit-title", HTML("&#10003; Crop-agnostic platform")),
            #   div(class = "benefit-text", "for diverse breeding programs")
            # ),
            # div(
            #   class = "benefit-card",
            #   div(class = "benefit-title", HTML("&#10003; Reproducible Research Workflows")),
            #   div(class = "benefit-text", "making all analysis repeatable")
            # )
          )
        ),

        # Footer Section
        div(
          class = "footer-section",
          # div(class = "footer-heading", "Project Funder"),
          # div(class = "footer-org", "Gates Foundation"),
          #
          # div(class = "footer-heading", "Lead Institution"),
          # div(class = "footer-org", "COLORADO STATE UNIVERSITY"),
          #
          # div(class = "logo-row",
          #     # Partner logos - uncomment and add paths when images are available
          #     # img(src = "logo_icar.png", class = "partner-logo"),
          #     # img(src = "logo_knust.png", class = "partner-logo"),
          #     # img(src = "logo_csir.png", class = "partner-logo"),
          #     # img(src = "logo_chibas.png", class = "partner-logo"),
          #     div(style = "color: #cbd5e0;", "[Partner Logos: www/logo_*.png]")
          # ),
          #
          # div(class = "logo-row",
          #     # img(src = "logo_cimmyt.png", class = "partner-logo"),
          #     # img(src = "logo_nebraska.png", class = "partner-logo"),
          #     # img(src = "logo_kansas.png", class = "partner-logo"),
          #     # img(src = "logo_hudsonalpha.png", class = "partner-logo"),
          #     div(style = "color: #cbd5e0;", "[Partner Logos Row 2]")
          # ),

          div(class = "copyright", HTML("&copy; 2025 panGenomeBreedr")),
          div(
            class = "developers",
            "Developed by: Alexander Wireko Kena, Israel Tawiah Tetteh, Cruet Burgos,
        Fanna Maina, Linly Banda, Jacques Faye, Benjamin Annor, Terry Felderhoff,
        Geoffrey Preston Norris"
          )
        )
      ),

      bslib::nav_item(), bslib::nav_item(), # space home tab from the rest of the tab

      ## Variant Annotation Tab
      bslib::nav_panel(
        title = "Variant Discovery",
        icon = icon("microscope"),
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
        title = "Marker Validation",
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
      )

      # # More tab for other information
      # navbarMenu(
      #   title = "More",
      #   icon = icon("ellipsis-v"),
      #
      #   # Tab for issues report
      #   bslib::nav_panel(
      #     title = "Report Issues",
      #   ),
      #   # Help tab
      #   bslib::nav_panel(title = "Help")
      # )
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
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
