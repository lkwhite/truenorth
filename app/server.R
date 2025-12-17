# server.R
# Compass Shiny app server logic

server <- function(input, output, session) {

  # ===========================================================================
  # Reactive values for app state
  # ===========================================================================

  values <- reactiveValues(
    trna_data = NULL,           # Loaded tRNA data frame with anticodon positions
    similarity_data = NULL,     # Pre-computed similarity matrix
    selected_ids = character(), # IDs of selected tRNAs
    selection_types = list(),   # Named list: ID -> "desired" or "avoid"
    probes = NULL               # Designed probes data frame
  )

  # ===========================================================================
  # Load data when organism changes
  # ===========================================================================

  observeEvent(input$organism, {
    # Show loading indicator
    showNotification("Loading tRNA data...", id = "loading", duration = NULL)

    # Load tRNA data
    values$trna_data <- load_trna_data(input$organism)

    # Load similarity data
    values$similarity_data <- load_sim_data(input$organism)

    # Clear selection when organism changes
    values$selected_ids <- character()
    values$selection_types <- list()
    values$probes <- NULL

    # Remove loading indicator
    removeNotification("loading")

    showNotification(
      paste0("Loaded ", nrow(values$trna_data), " tRNAs for ",
             names(ORGANISMS)[ORGANISMS == input$organism]),
      type = "message",
      duration = 3
    )
  }, ignoreNULL = FALSE)

  # ===========================================================================
  # Update selection summary in header
  # ===========================================================================

  observe({
    n_desired <- sum(unlist(values$selection_types) == "desired")
    n_avoid <- sum(unlist(values$selection_types) == "avoid")

    # Update the summary display using JavaScript
    session$sendCustomMessage("updateSelectionSummary", list(
      desired = n_desired,
      avoid = n_avoid
    ))
  })

  # Add JavaScript handler for updating selection summary
  tags$script(HTML("
    Shiny.addCustomMessageHandler('updateSelectionSummary', function(data) {
      document.getElementById('desired-count').innerText = data.desired + ' desired';
      document.getElementById('avoid-count').innerText = data.avoid + ' avoid';
    });
  "))

  # ===========================================================================
  # Module servers
  # ===========================================================================

  # Browse module - returns updated selection
  browse_result <- browseServer("browse", values)

  # Update values when browse module changes selection
  observeEvent(browse_result$selected_ids(), {
    values$selected_ids <- browse_result$selected_ids()
  }, ignoreNULL = FALSE)

  observeEvent(browse_result$selection_types(), {
    values$selection_types <- browse_result$selection_types()
  }, ignoreNULL = FALSE)

  # Design module - returns probes
  design_result <- designServer("design", values)

  # Update probes when design completes
  observeEvent(design_result$probes(), {
    values$probes <- design_result$probes()
    if (!is.null(values$probes) && nrow(values$probes) > 0) {
      # Navigate to results tab
      updateNavbarPage(session, "navbar", selected = "Results")
    }
  }, ignoreNULL = FALSE)

  # Results module
  resultsServer("results", values)

}
