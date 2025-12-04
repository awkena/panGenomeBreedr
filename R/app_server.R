#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#'
app_server <- function(input, output, session) {

  # Auto-stop R when the browser closes
  session$onSessionEnded(function(){
    stopApp()
  })

  #------ server side for site redirecting-------#
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
  observeEvent(input$refresh_btn, {
    # Show confirmation modal
    showModal(
      modalDialog(
        title = "Refresh Application",
        footer = NULL, # Remove default buttons
        tagList(
          p("This will refresh the entire application and reset all data."),
          p("Do you wish to continue?"),

          #  Choice
          radioButtons(
            inputId = "refresh_choice",
            label = NULL,
            choices = c("Yes, refresh the application" = TRUE),
            selected = character(0)
          ),
          div(
            style = "text-align: right;",
            actionButton("confirm_refresh", "Confirm",
                         class = "btn-secondary"
            ),
            actionButton("cancel_refresh", "Cancel",
                         class = "btn-danger"
            )
          )
        )
      )
    )
  })

  # Handle confirm button
  observeEvent(input$confirm_refresh, {
    if (!is.null(input$refresh_choice) && input$refresh_choice == TRUE) {
      removeModal()
      shinyWidgets::show_toast(
        title = "Refreshing application...",
        type = "info",
        timer = 300
      )
      session$reload() # refresh entire web app
    }
  })

  # Remove modal when cancel is clicked
  observeEvent(input$cancel_refresh, {
    removeModal()
  })

  #-------

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
    #--Knitted with qc plot
    # # server side for  plate layout
    # mod_mv_plate_plot_server("mv_plate_plot_1",
    #                          kasp_data = import_data_entities(),
    #                          color_coded = color_code_res()
    # )
  })

  ## -------------------- Decision support, server-side ------------------------##
  mod_ds_trait_introg_hypothesis_test_server("ds_trait_introg_hypothesis_test_1")

  mod_ds_marker_ass_bac_server("ds_marker_ass_bac_1") # marker assisted back-cross

  mod_ds_foreground_select_server("ds_foreground_select_1") # foreground
}
