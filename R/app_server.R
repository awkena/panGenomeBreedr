#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#'
app_server <- function(input, output, session) {

  #------ server side for site redirecting to vignette-------#
  #- variant discovery
  observeEvent(input$btn_variant_discovery,{
  session$sendCustomMessage('open_link',
                            list(url = 'https://awkena.github.io/panGenomeBreedr/articles/panGenomeBreedr_Workflows.html#variant-discovery'))
  })

  #- kasp marker design.
  observeEvent(input$btn_kasp_marker,{
    session$sendCustomMessage('open_link',
                              list(url = 'https://awkena.github.io/panGenomeBreedr/articles/panGenomeBreedr_Workflows.html#kasp-marker-design'))
  })

  #- marker validation
  observeEvent(input$btn_marker_validation,{
    session$sendCustomMessage('open_link',
                              list(url = 'https://awkena.github.io/panGenomeBreedr/articles/panGenomeBreedr_Workflows.html#kasp-marker-validation'))
  })

  #- decision support tool
  observeEvent(input$btn_introgression,{
    session$sendCustomMessage('open_link',
                              list(url = 'https://awkena.github.io/panGenomeBreedr/articles/panGenomeBreedr_Workflows.html#decision-support-for-trait-introgression-and-mabc'))
  })

  #-------- Let user referesh application.
  observeEvent(input$power_btn, {
    showModal(
      modalDialog(
        title = div(
          style = "display: flex; align-items: center;",
          icon("power-off", class = "text-danger me-2"),
          tags$h4("System Controls", class = "mb-0")
        ),
        footer = NULL,
        tagList(
          p("What would you like to do with the application?"),

          div(
            class = "d-grid gap-3", # Bootstrap vertical stack

            # OPTION 1: REFRESH
            actionButton(
              "btn_modal_refresh",
              tagList(icon("sync"), "Refresh Application"),
              class = "btn-outline-primary p-3 text-start",
              style = "border-radius: 10px;"
            ),
            p(class = "text-muted small ms-2", "Resets all data and reloads the current page."),

            # OPTION 2: SHUT DOWN
            actionButton(
              "btn_modal_shutdown",
              tagList(icon("stop-circle"), "Shut Down Application"),
              class = "btn-outline-danger p-3 text-start",
              style = "border-radius: 10px;"
            ),
            p(class = "text-muted small ms-2", "Closes the app and stops the R process.")
          ),

          hr(),
          div(class = "text-end",
              actionButton("btn_modal_cancel", "Cancel", class = "btn-secondary btn-sm")
          )
        )
      )
    )
  })

  # --- Execution Logic for refresh & shutdown ---
  # Handle Refresh
  observeEvent(input$btn_modal_refresh, {
    removeModal()
    session$reload()
  })

  # Handle Shutdown
  observeEvent(input$btn_modal_shutdown, {
    removeModal()
    stopApp()
  })

  # Handle Cancel
  observeEvent(input$btn_modal_cancel, {
    removeModal()
  })



  #--------- Variant discovery server side --------------#
  mod_variant_discovery_server("variant_discovery_1")


  #-------------- kasp marker design server side-----------#
  mod_kasp_marker_design_server(id = "kasp_marker_design_1")


  ## ------------------------- Marker Validation Server side -------------------##
  # import_data_entities; returns kasp data after being read
  import_data_entities <- mod_mv_read_kasp_csv_server("mv_read_kasp_csv_1")


  # nsamples & alleles
  observe({
    req(import_data_entities())
    mod_mv_nsamples_plate_server("mv_nsamples_plate_1",
                                 kasp_data = import_data_entities()
    )
    mod_mv_get_alleles_server("mv_get_alleles_1",
                              kasp_data = import_data_entities()
    )
  })

  # Color coding server side
  color_code_res <- reactiveVal(NULL) # Create a reactive null to hold result.

  observe({
    req(import_data_entities())
    color_code_res(
      mod_mv_kasp_color_server("mv_kasp_color_1",
                               kasp_data = import_data_entities()
      )
    )
  })

  # Predictive summary and status server side
  observe({
    req(import_data_entities(), color_code_res())

    mod_mv_pred_sum_stat_server("mv_pred_sum_stat_1",
                                color_code_res = color_code_res(),
                                kasp_data = import_data_entities()
    )
  })

  # Qc plot server portion
  observe({
    req(color_code_res(), import_data_entities())
    # Server side for qcplot
    mod_mv_kasp_qc_ggplot_server(
      id = "mv_kasp_qc_ggplot_1",
      kasp_data = import_data_entities(),
      color_coded = color_code_res()
    )
  })

  ## -------------------- Decision support, server-side ------------------------##
  mod_ds_trait_introg_hypothesis_test_server("ds_trait_introg_hypothesis_test_1")

  mod_ds_marker_ass_bac_server("ds_marker_ass_bac_1") # marker assisted back-cross

  mod_ds_foreground_select_server("ds_foreground_select_1") # foreground
}
