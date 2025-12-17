# wizard_step4.R
# Step 4: Design Results - probe table and export

wizardStep4UI <- function(id) {
  ns <- NS(id)

  tagList(
    # Selection context
    uiOutput(ns("selection_context")),

    tags$div(
      class = "text-center mb-4",
      tags$h3("Designed Probes"),
      tags$p(class = "text-muted",
             "Review your probes and export for synthesis")
    ),

    uiOutput(ns("results_content"))
  )
}

wizardStep4Server <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Selection context banner
    output$selection_context <- renderUI({
      req(values$wizard_goal)
      req(values$wizard_selection$ids)

      n_targets <- length(values$wizard_selection$ids)
      goal_label <- switch(values$wizard_goal,
        "specific" = "Specific tRNAs",
        "isoacceptor" = paste0("Isoacceptor: ", values$wizard_selection$anticodon),
        "amino_acid" = paste0("Amino acid: ", values$wizard_selection$amino_acid),
        "consensus" = "All tRNAs (consensus)"
      )

      tags$div(
        class = "alert alert-light border mb-4",
        tags$div(
          class = "d-flex justify-content-between align-items-center",
          tags$div(
            tags$strong("Target: "), goal_label,
            tags$span(class = "badge bg-primary ms-2", paste(n_targets, "tRNAs"))
          ),
          actionLink(ns("change_selection"), "Change selection", class = "small")
        )
      )
    })

    # Handle change selection link
    observeEvent(input$change_selection, {
      values$wizard_step <- 2
    })

    output$results_content <- renderUI({
      probes <- values$wizard_probes

      if (is.null(probes)) {
        return(tags$div(
          class = "text-center text-muted",
          tags$p("Designing probes..."),
          tags$div(class = "spinner-border", role = "status")
        ))
      }

      if (nrow(probes) == 0) {
        return(tags$div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          tags$strong(" No probes found"),
          tags$p("Could not design probes meeting the specified criteria. ",
                 "Try adjusting the parameters or selecting different targets.")
        ))
      }

      # Summary stats
      n_probes <- nrow(probes)
      avg_tm <- mean(probes$tm_nn, na.rm = TRUE)
      avg_gc <- mean(probes$gc_content, na.rm = TRUE)

      tagList(
        # Summary cards
        fluidRow(
          column(
            width = 3,
            tags$div(
              class = "card text-center",
              tags$div(
                class = "card-body",
                tags$h2(style = "color: #0072B2;", n_probes),
                tags$small(class = "text-muted", "Probes designed")
              )
            )
          ),
          column(
            width = 3,
            tags$div(
              class = "card text-center",
              tags$div(
                class = "card-body",
                tags$h2(sprintf("%.1f\u00B0C", avg_tm)),
                tags$small(class = "text-muted", "Average Tm")
              )
            )
          ),
          column(
            width = 3,
            tags$div(
              class = "card text-center",
              tags$div(
                class = "card-body",
                tags$h2(sprintf("%.0f%%", avg_gc)),
                tags$small(class = "text-muted", "Average GC")
              )
            )
          ),
          column(
            width = 3,
            tags$div(
              class = "card text-center",
              tags$div(
                class = "card-body",
                downloadButton(ns("download_probes"), "Export CSV",
                               class = "btn-primary")
              )
            )
          )
        ),

        tags$hr(),

        # Probe table
        tags$div(
          class = "card",
          tags$div(
            class = "card-header bg-light d-flex justify-content-between align-items-center",
            tags$strong("Probe Details"),
            tags$div(
              class = "btn-group btn-group-sm",
              actionButton(ns("show_all"), "All", class = "btn-outline-secondary"),
              actionButton(ns("show_best"), "Best Only", class = "btn-outline-secondary active")
            )
          ),
          tags$div(
            class = "card-body",
            DT::dataTableOutput(ns("probe_table"))
          )
        ),

        # Probe detail view
        tags$div(
          class = "card mt-4",
          tags$div(
            class = "card-header bg-light",
            tags$strong("Probe Detail")
          ),
          tags$div(
            class = "card-body",
            uiOutput(ns("probe_detail"))
          )
        )
      )
    })

    # Probe table
    output$probe_table <- DT::renderDataTable({
      req(values$wizard_probes)
      probes <- values$wizard_probes

      # Format for display using actual column names
      display_df <- data.frame(
        Rank = probes$rank,
        Sequence = probes$probe_sequence,
        Region = probes$trna_region,
        Position = paste0(probes$start, "-", probes$end),
        Tm = sprintf("%.1f", probes$tm_nn),
        GC = sprintf("%.0f%%", probes$gc_content),
        Quality = probes$quality,
        stringsAsFactors = FALSE
      )

      # Add warning icon for anticodon overlap
      if ("overlaps_anticodon" %in% names(probes)) {
        display_df$Region <- ifelse(
          probes$overlaps_anticodon,
          paste0(probes$trna_region, " \u26A0"),
          probes$trna_region
        )
      }

      DT::datatable(
        display_df,
        selection = "single",
        options = list(
          pageLength = 10,
          dom = "tip",
          columnDefs = list(
            list(className = "dt-center", targets = c(3, 4, 5, 6))
          )
        ),
        rownames = FALSE,
        escape = FALSE  # Allow HTML/unicode in cells
      )
    })

    # Probe detail view
    output$probe_detail <- renderUI({
      req(values$wizard_probes)
      selected <- input$probe_table_rows_selected

      if (is.null(selected) || length(selected) == 0) {
        return(tags$p(class = "text-muted text-center",
                      "Click on a probe in the table above to see details"))
      }

      probe <- values$wizard_probes[selected, ]

      # Check for anticodon overlap
      has_anticodon_overlap <- if ("overlaps_anticodon" %in% names(probe)) {
        probe$overlaps_anticodon
      } else {
        FALSE
      }

      tagList(
        # Anticodon warning if applicable
        if (has_anticodon_overlap) {
          tags$div(
            class = "alert alert-warning py-2 mb-3",
            icon("exclamation-triangle"),
            tags$strong(" Note: "),
            "This probe overlaps the anticodon region, which is heavily modified. ",
            "Hybridization efficiency may be reduced."
          )
        },

        fluidRow(
          column(
            width = 6,
            tags$h5("Probe Sequence"),
            tags$div(
              class = "p-3 bg-light rounded",
              style = "font-family: monospace; font-size: 1.2em; word-break: break-all;",
              probe$probe_sequence
            ),
            tags$div(
              class = "mt-2",
              tags$small(class = "text-muted",
                         "5' \u2192 3' orientation, ready for synthesis")
            )
          ),
          column(
            width = 6,
            tags$h5("Properties"),
            tags$table(
              class = "table table-sm",
              tags$tbody(
                tags$tr(
                  tags$td("tRNA Region:"),
                  tags$td(tags$strong(probe$trna_region))
                ),
                tags$tr(
                  tags$td("Position:"),
                  tags$td(tags$strong(paste0(probe$start, "-", probe$end)))
                ),
                tags$tr(
                  tags$td("Length:"),
                  tags$td(tags$strong(probe$length, " nt"))
                ),
                tags$tr(
                  tags$td("Melting Temp:"),
                  tags$td(tags$strong(sprintf("%.1f\u00B0C", probe$tm_nn)))
                ),
                tags$tr(
                  tags$td("GC Content:"),
                  tags$td(tags$strong(sprintf("%.1f%%", probe$gc_content)))
                ),
                tags$tr(
                  tags$td("Quality:"),
                  tags$td(tags$strong(probe$quality))
                ),
                tags$tr(
                  tags$td("Target Conservation:"),
                  tags$td(tags$strong(sprintf("%.1f%%", probe$desired_conservation)))
                )
              )
            )
          )
        ),

        # Reference sequence info
        tags$div(
          class = "mt-3",
          tags$small(class = "text-muted",
                     "Reference: ", format_trna_id(probe$reference_id))
        ),

        # Hybridization diagram - Targets
        tags$hr(),
        tags$div(
          class = "d-flex justify-content-between align-items-center mb-2",
          tags$h5(class = "mb-0", "Target Binding"),
          tags$div(
            class = "form-check form-switch",
            tags$input(
              type = "checkbox",
              class = "form-check-input",
              id = ns("show_all_targets"),
              onclick = sprintf("Shiny.setInputValue('%s', this.checked, {priority: 'event'})", ns("show_all_targets"))
            ),
            tags$label(class = "form-check-label small", `for` = ns("show_all_targets"), "Show all")
          )
        ),
        uiOutput(ns("target_binding_diagram")),

        # Off-target analysis (check ALL compartments for cross-hybridization)
        tags$div(
          class = "d-flex justify-content-between align-items-center mb-2 mt-4",
          tags$h5(class = "mb-0", "Potential Off-Targets"),
          tags$div(
            class = "form-check form-switch",
            tags$input(
              type = "checkbox",
              class = "form-check-input",
              id = ns("show_all_offtargets"),
              onclick = sprintf("Shiny.setInputValue('%s', this.checked, {priority: 'event'})", ns("show_all_offtargets"))
            ),
            tags$label(class = "form-check-label small", `for` = ns("show_all_offtargets"), "Show all")
          )
        ),
        uiOutput(ns("offtarget_diagram"))
      )
    })

    # Reactive target binding diagram with toggle
    output$target_binding_diagram <- renderUI({
      req(values$wizard_probes)

      selected_row <- input$probe_table_rows_selected
      if (is.null(selected_row) || length(selected_row) == 0) {
        return(NULL)
      }

      probe <- values$wizard_probes[selected_row, ]
      target_ids <- values$wizard_selection$ids
      targets_df <- values$trna_data[values$trna_data$id %in% target_ids, ]

      if (nrow(targets_df) == 0) {
        return(tags$div(class = "text-muted", "Target sequence not available"))
      }

      # Use toggle to determine max_show (NULL = show all)
      show_all <- input$show_all_targets %||% FALSE
      max_show <- if (show_all) Inf else 5

      render_multi_target_hybridization(
        targets_df = targets_df,
        probe_sequence = probe$probe_sequence,
        start = probe$start,
        end = probe$end,
        reference_id = probe$reference_id,
        max_show = max_show
      )
    })

    # Reactive off-target diagram with toggle
    output$offtarget_diagram <- renderUI({
      req(values$wizard_probes)

      selected_row <- input$probe_table_rows_selected
      if (is.null(selected_row) || length(selected_row) == 0) {
        return(NULL)
      }

      probe <- values$wizard_probes[selected_row, ]
      target_ids <- values$wizard_selection$ids
      all_trnas <- if (!is.null(values$trna_data_all)) values$trna_data_all else values$trna_data
      off_target_df <- all_trnas[!all_trnas$id %in% target_ids, ]

      if (nrow(off_target_df) == 0) {
        return(tags$div(
          class = "alert alert-success py-2",
          icon("check-circle"),
          " No off-targets - all tRNAs are targets"
        ))
      }

      # Calculate mismatches for each off-target
      off_target_df$binding_region <- substr(off_target_df$sequence, probe$start, probe$end)
      probe_reversed <- paste(rev(strsplit(probe$probe_sequence, "")[[1]]), collapse = "")
      probe_chars <- strsplit(probe_reversed, "")[[1]]

      off_target_df$n_mismatches <- sapply(off_target_df$binding_region, function(region) {
        region_chars <- strsplit(region, "")[[1]]
        mismatches <- 0
        for (i in seq_along(region_chars)) {
          target_base <- toupper(region_chars[i])
          probe_base <- if (i <= length(probe_chars)) toupper(probe_chars[i]) else "?"
          is_match <- (target_base == "A" && probe_base == "T") ||
                      (target_base == "T" && probe_base == "A") ||
                      (target_base == "G" && probe_base == "C") ||
                      (target_base == "C" && probe_base == "G")
          if (!is_match) mismatches <- mismatches + 1
        }
        mismatches
      })

      # Sort by mismatches
      off_target_df <- off_target_df[order(off_target_df$n_mismatches), ]

      # Use toggle to determine max_show
      show_all <- input$show_all_offtargets %||% FALSE
      max_show <- if (show_all) Inf else 3
      closest_off_targets <- if (show_all) off_target_df else head(off_target_df, max_show)

      min_mismatches <- min(closest_off_targets$n_mismatches)

      # Specificity alert
      alert <- if (min_mismatches >= 5) {
        tags$div(
          class = "alert alert-success py-2 mb-3",
          icon("check-circle"),
          paste0(" Good specificity: closest off-target has ", min_mismatches, " mismatches")
        )
      } else if (min_mismatches >= 3) {
        tags$div(
          class = "alert alert-warning py-2 mb-3",
          icon("exclamation-triangle"),
          paste0(" Moderate specificity: closest off-target has ", min_mismatches, " mismatches")
        )
      } else {
        tags$div(
          class = "alert alert-danger py-2 mb-3",
          icon("exclamation-circle"),
          paste0(" Low specificity: closest off-target has only ", min_mismatches,
                 " mismatch", if (min_mismatches != 1) "es" else "")
        )
      }

      tagList(
        alert,
        render_multi_target_hybridization(
          targets_df = closest_off_targets,
          probe_sequence = probe$probe_sequence,
          start = probe$start,
          end = probe$end,
          reference_id = NULL,
          max_show = max_show
        )
      )
    })

    # Download handler
    output$download_probes <- downloadHandler(
      filename = function() {
        goal <- values$wizard_goal %||% "probes"
        target <- switch(goal,
          "amino_acid" = values$wizard_selection$amino_acid,
          "isoacceptor" = values$wizard_selection$anticodon,
          "specific" = "selected",
          "consensus" = "pantRNA",
          "probes"
        )
        paste0("compass_", target, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        probes <- values$wizard_probes

        # Format for export with actual column names
        export_df <- data.frame(
          rank = probes$rank,
          probe_sequence = probes$probe_sequence,
          trna_region = probes$trna_region,
          start = probes$start,
          end = probes$end,
          length = probes$length,
          tm_celsius = round(probes$tm_nn, 1),
          gc_percent = round(probes$gc_content, 1),
          quality = probes$quality,
          target_conservation = round(probes$desired_conservation, 1),
          selectivity_score = round(probes$selectivity_score, 1),
          overlaps_anticodon = probes$overlaps_anticodon,
          reference_id = probes$reference_id,
          stringsAsFactors = FALSE
        )

        write.csv(export_df, file, row.names = FALSE)
      }
    )
  })
}
