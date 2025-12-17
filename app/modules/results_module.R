# results_module.R
# Shiny module for displaying probe design results

# =============================================================================
# UI
# =============================================================================

resultsUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        # Header with export button
        tags$div(
          class = "d-flex justify-content-between align-items-center mb-3",
          tags$h4("Designed Probes"),
          downloadButton(ns("export_csv"), "Export CSV", class = "btn-outline-primary")
        ),

        # Results or placeholder
        uiOutput(ns("results_content"))
      )
    )
  )
}

# =============================================================================
# Server
# =============================================================================

resultsServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # -------------------------------------------------------------------------
    # Results content
    # -------------------------------------------------------------------------
    output$results_content <- renderUI({
      if (is.null(values$probes) || nrow(values$probes) == 0) {
        return(tags$div(
          class = "alert alert-info",
          icon("info-circle"),
          " No probes designed yet. Go to ",
          tags$strong("Design Probes"),
          " to create probes for your selected targets."
        ))
      }

      # Show the table
      tagList(
        # Summary
        tags$div(
          class = "alert alert-success",
          icon("check-circle"),
          paste0(" ", nrow(values$probes), " probes designed. ",
                 "Click on a row for details.")
        ),

        # Table
        DTOutput(ns("probe_table"))
      )
    })

    # -------------------------------------------------------------------------
    # Probe table
    # -------------------------------------------------------------------------
    output$probe_table <- renderDT({
      req(values$probes)

      # Select and format columns for display
      df <- values$probes

      # Determine which columns exist
      display_cols <- c("rank", "probe_sequence", "start", "end", "length",
                        "gc_content", "tm_nn")

      # Add optional columns if they exist
      if ("desired_conservation" %in% names(df)) {
        display_cols <- c(display_cols, "desired_conservation")
      }
      if ("avoid_divergence" %in% names(df)) {
        display_cols <- c(display_cols, "avoid_divergence")
      }
      if ("quality" %in% names(df)) {
        display_cols <- c(display_cols, "quality")
      }

      # Filter to existing columns
      display_cols <- display_cols[display_cols %in% names(df)]
      df <- df[, display_cols, drop = FALSE]

      # Rename for display
      col_names <- c(
        rank = "Rank",
        probe_sequence = "Probe Sequence (5'→3')",
        start = "Start",
        end = "End",
        length = "Length",
        gc_content = "GC%",
        tm_nn = "Tm (°C)",
        desired_conservation = "Conservation",
        avoid_divergence = "Divergence",
        quality = "Quality"
      )

      names(df) <- col_names[names(df)]

      # Format numeric columns
      if ("GC%" %in% names(df)) {
        df$`GC%` <- round(df$`GC%`, 1)
      }
      if ("Tm (°C)" %in% names(df)) {
        df$`Tm (°C)` <- round(df$`Tm (°C)`, 1)
      }
      if ("Conservation" %in% names(df)) {
        df$Conservation <- round(df$Conservation, 1)
      }
      if ("Divergence" %in% names(df)) {
        df$Divergence <- round(df$Divergence, 1)
      }

      datatable(
        df,
        options = list(
          pageLength = 10,
          dom = 'frtip',
          scrollX = TRUE,
          columnDefs = list(
            list(className = 'dt-center', targets = '_all')
          )
        ),
        selection = 'single',
        rownames = FALSE,
        class = 'display compact'
      ) %>%
        formatStyle(
          'Probe Sequence (5\'→3\')',
          fontFamily = 'monospace',
          fontSize = '0.9em'
        )
    })

    # -------------------------------------------------------------------------
    # Export CSV
    # -------------------------------------------------------------------------
    output$export_csv <- downloadHandler(
      filename = function() {
        paste0("compass_probes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        req(values$probes)
        write.csv(values$probes, file, row.names = FALSE)
      }
    )

  })
}
