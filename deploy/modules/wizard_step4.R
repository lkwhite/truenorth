# wizard_step4.R
# Step 4: Design Results - probe table with coverage selection and export

wizardStep4UI <- function(id) {
  ns <- NS(id)

  tagList(
    # Selection context
    uiOutput(ns("selection_context")),

    tags$div(
      class = "text-center mb-4",
      tags$h3("Designed Probes"),
      tags$p(class = "text-muted",
             "Select probes to build your probe set and track target coverage")
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
        # Summary cards row
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
                class = "card-body py-2",
                tags$small(class = "text-muted d-block mb-1", "Export selected"),
                downloadButton(ns("download_probes"), "Export CSV",
                               class = "btn-primary btn-sm")
              )
            )
          )
        ),

        # Coverage tracker card
        tags$div(
          class = "card mt-3 mb-3 border-primary",
          tags$div(
            class = "card-header bg-primary text-white d-flex justify-content-between align-items-center",
            tags$span(icon("vials"), " Your Probe Set"),
            tags$span(class = "badge bg-light text-primary", uiOutput(ns("coverage_badge"), inline = TRUE))
          ),
          tags$div(
            class = "card-body",
            uiOutput(ns("coverage_tracker"))
          )
        ),

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

    # Probe table with multi-select (click to select for probe set AND view details)
    output$probe_table <- DT::renderDataTable({
      req(values$wizard_probes)
      probes <- values$wizard_probes
      coverage_data <- values$coverage_data

      # Format coverage columns
      if (!is.null(coverage_data) && "n_targets_hit" %in% names(probes)) {
        n_targets <- coverage_data$n_targets
        covers_col <- sprintf("%d/%d (%.0f%%)", probes$n_targets_hit, n_targets, probes$coverage_pct)
      } else {
        covers_col <- rep("-", nrow(probes))
      }

      # Format for display
      display_df <- data.frame(
        Rank = probes$rank,
        Sequence = probes$probe_sequence,
        Region = probes$trna_region,
        Covers = covers_col,
        Tm = sprintf("%.1f", probes$tm_nn),
        Quality = probes$quality,
        stringsAsFactors = FALSE
      )

      # Add warning icon for modification risk
      if ("modification_penalty" %in% names(probes)) {
        display_df$Region <- ifelse(
          probes$modification_penalty >= 20,
          paste0(probes$trna_region, " \u26A0\u26A0"),  # High risk: double warning
          ifelse(
            probes$modification_penalty >= 10,
            paste0(probes$trna_region, " \u26A0"),  # Moderate risk: single warning
            probes$trna_region
          )
        )
      } else if ("overlaps_anticodon" %in% names(probes)) {
        # Legacy fallback
        display_df$Region <- ifelse(
          probes$overlaps_anticodon,
          paste0(probes$trna_region, " \u26A0"),
          probes$trna_region
        )
      }

      DT::datatable(
        display_df,
        selection = list(mode = "multiple", selected = which(probes$rank %in% values$selected_probes)),
        options = list(
          pageLength = 10,
          dom = "tip",
          columnDefs = list(
            list(className = "dt-center", targets = c(0, 3, 4, 5))
          )
        ),
        rownames = FALSE,
        escape = FALSE  # Allow HTML/unicode in cells
      )
    })

    # Sync table selection with values$selected_probes
    observeEvent(input$probe_table_rows_selected, {
      req(values$wizard_probes)
      selected_rows <- input$probe_table_rows_selected
      if (is.null(selected_rows)) {
        values$selected_probes <- integer()
      } else {
        values$selected_probes <- values$wizard_probes$rank[selected_rows]
      }
    }, ignoreNULL = FALSE)

    # Coverage badge (shows in card header)
    output$coverage_badge <- renderUI({
      coverage_data <- values$coverage_data
      selected <- values$selected_probes %||% integer()

      if (is.null(coverage_data) || length(selected) == 0) {
        return(tags$span("0 probes"))
      }

      coverage <- get_cumulative_coverage(selected, coverage_data$matrix, coverage_data$target_ids)
      tags$span(sprintf("%d probe%s", length(selected), if (length(selected) == 1) "" else "s"))
    })

    # Coverage tracker content
    output$coverage_tracker <- renderUI({
      coverage_data <- values$coverage_data
      selected <- values$selected_probes %||% integer()

      if (is.null(coverage_data)) {
        return(tags$p(class = "text-muted mb-0", "Coverage data not available"))
      }

      n_targets <- coverage_data$n_targets
      target_ids <- coverage_data$target_ids

      if (length(selected) == 0) {
        # No probes selected - show recommendation
        estimate <- values$coverage_estimate
        if (!is.null(estimate) && nrow(estimate) > 0) {
          full_coverage_row <- estimate[estimate$coverage_pct >= 99.9, ]
          if (nrow(full_coverage_row) > 0) {
            n_for_full <- min(full_coverage_row$n_probes)
            rec_ranks <- full_coverage_row$probe_ranks[full_coverage_row$n_probes == n_for_full][[1]]
            return(tagList(
              tags$p(class = "text-muted mb-2",
                     icon("lightbulb"), " Select probes below to build your probe set"),
              tags$div(
                class = "alert alert-info py-2 mb-0",
                tags$strong("Recommendation: "),
                sprintf("Select %d probe%s for full coverage (ranks: %s)",
                        n_for_full,
                        if (n_for_full == 1) "" else "s",
                        paste(rec_ranks, collapse = ", ")),
                tags$br(),
                actionLink(ns("select_recommended"), "Select recommended probes",
                           class = "small")
              )
            ))
          }
        }
        return(tags$p(class = "text-muted mb-0",
                      icon("lightbulb"), " Select probes below to build your probe set"))
      }

      # Get cumulative coverage for selected probes
      coverage <- get_cumulative_coverage(selected, coverage_data$matrix, target_ids)
      pct <- coverage$coverage_pct

      # Progress bar color based on coverage
      bar_class <- if (pct >= 100) "bg-success" else if (pct >= 80) "bg-info" else if (pct >= 50) "bg-warning" else "bg-danger"

      tagList(
        # Coverage progress bar
        tags$div(
          class = "mb-3",
          tags$div(
            class = "d-flex justify-content-between mb-1",
            tags$span("Target Coverage"),
            tags$strong(sprintf("%d/%d (%.0f%%)", coverage$n_covered, n_targets, pct))
          ),
          tags$div(
            class = "progress",
            style = "height: 20px;",
            tags$div(
              class = paste("progress-bar", bar_class),
              role = "progressbar",
              style = sprintf("width: %.0f%%;", pct),
              `aria-valuenow` = round(pct),
              `aria-valuemin` = 0,
              `aria-valuemax` = 100
            )
          )
        ),

        # Covered/uncovered lists
        fluidRow(
          column(
            width = 6,
            tags$div(
              class = "small",
              tags$strong(class = "text-success", icon("check-circle"), " Covered: "),
              if (coverage$n_covered > 0) {
                tags$span(
                  title = paste(format_trna_id(coverage$covered_ids), collapse = ", "),
                  if (coverage$n_covered <= 5) {
                    paste(format_trna_id(coverage$covered_ids), collapse = ", ")
                  } else {
                    paste0(paste(format_trna_id(head(coverage$covered_ids, 4)), collapse = ", "),
                           sprintf(" +%d more", coverage$n_covered - 4))
                  }
                )
              } else {
                tags$span(class = "text-muted", "None")
              }
            )
          ),
          column(
            width = 6,
            tags$div(
              class = "small",
              if (coverage$n_uncovered > 0) {
                tagList(
                  tags$strong(class = "text-danger", icon("times-circle"), " Not covered: "),
                  tags$span(
                    title = paste(format_trna_id(coverage$uncovered_ids), collapse = ", "),
                    if (coverage$n_uncovered <= 5) {
                      paste(format_trna_id(coverage$uncovered_ids), collapse = ", ")
                    } else {
                      paste0(paste(format_trna_id(head(coverage$uncovered_ids, 4)), collapse = ", "),
                             sprintf(" +%d more", coverage$n_uncovered - 4))
                    }
                  )
                )
              } else {
                tags$span(class = "text-success", icon("check-circle"), " All targets covered!")
              }
            )
          )
        ),

        # Clear selection link
        if (length(selected) > 0) {
          tags$div(
            class = "mt-2 text-end",
            actionLink(ns("clear_selection"), "Clear selection", class = "small text-muted")
          )
        }
      )
    })

    # Handle "Select recommended probes" link
    observeEvent(input$select_recommended, {
      estimate <- values$coverage_estimate
      if (!is.null(estimate) && nrow(estimate) > 0) {
        full_coverage_row <- estimate[estimate$coverage_pct >= 99.9, ]
        if (nrow(full_coverage_row) > 0) {
          n_for_full <- min(full_coverage_row$n_probes)
          rec_ranks <- full_coverage_row$probe_ranks[full_coverage_row$n_probes == n_for_full][[1]]
          values$selected_probes <- rec_ranks
        }
      }
    })

    # Handle clear selection link
    observeEvent(input$clear_selection, {
      values$selected_probes <- integer()
    })

    # Probe detail view - shows details for last clicked probe
    output$probe_detail <- renderUI({
      req(values$wizard_probes)

      # Use last clicked row for details (works with multi-select)
      last_clicked <- input$probe_table_row_last_clicked
      selected <- input$probe_table_rows_selected

      if (is.null(last_clicked) && (is.null(selected) || length(selected) == 0)) {
        return(tags$p(class = "text-muted text-center",
                      "Click on a probe to select it and see details. ",
                      "Ctrl+click or Shift+click to select multiple probes."))
      }

      # Show details for last clicked, or first selected if no click tracked
      detail_row <- if (!is.null(last_clicked)) last_clicked else selected[1]
      probe <- values$wizard_probes[detail_row, ]

      # Check modification penalty level
      mod_penalty <- if ("modification_penalty" %in% names(probe)) {
        probe$modification_penalty
      } else if ("overlaps_anticodon" %in% names(probe) && probe$overlaps_anticodon) {
        25  # Legacy fallback
      } else {
        0
      }

      # Determine modification risk level
      mod_risk <- if (mod_penalty >= 20) {
        list(level = "high", class = "alert-warning",
             message = "This probe overlaps heavily modified regions (anticodon loop). Hybridization efficiency may be significantly reduced.")
      } else if (mod_penalty >= 10) {
        list(level = "moderate", class = "alert-info",
             message = "This probe overlaps moderately modified regions (D-loop or T\u03A8C loop). Some impact on hybridization possible.")
      } else {
        NULL
      }

      tagList(
        # Modification warning if applicable
        if (!is.null(mod_risk)) {
          tags$div(
            class = paste("alert py-2 mb-3", mod_risk$class),
            icon("exclamation-triangle"),
            tags$strong(" Modification note: "),
            mod_risk$message
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
                  tags$td(
                    tags$strong(sprintf("%.1f\u00B0C", probe$tm_nn)),
                    tags$sup(
                      style = "color: #E69F00; cursor: help;",
                      title = "Calculated for unmodified DNA/RNA. tRNA modifications may reduce effective Tm by 5-15\u00B0C.",
                      " *"
                    )
                  )
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

        # Reference sequence info and Tm footnote
        tags$div(
          class = "mt-3",
          tags$small(class = "text-muted",
                     "Reference: ", format_trna_id(probe$reference_id)),
          tags$br(),
          tags$small(
            style = "color: #E69F00;",
            "* Tm calculated for unmodified sequences. tRNA bases are heavily modified, ",
            "especially in the anticodon loop. Actual hybridization efficiency may vary."
          )
        ),

        # Hybridization diagram - Targets
        tags$hr(),
        tags$div(
          class = "d-flex justify-content-between align-items-center mb-2",
          tags$h5(class = "mb-0", "Target Binding"),
          tags$div(
            class = "d-flex align-items-center gap-2",
            tags$span(class = "small text-muted", "Top 5"),
            tags$div(
              class = "form-check form-switch mb-0",
              tags$input(
                type = "checkbox",
                class = "form-check-input",
                id = ns("show_all_targets"),
                onclick = sprintf("Shiny.setInputValue('%s', this.checked, {priority: 'event'})", ns("show_all_targets"))
              ),
              tags$label(class = "form-check-label small", `for` = ns("show_all_targets"), "Show all")
            )
          )
        ),
        uiOutput(ns("target_binding_diagram")),

        # Off-target analysis (check ALL compartments for cross-hybridization)
        tags$div(
          class = "d-flex justify-content-between align-items-center mb-2 mt-4",
          tags$h5(class = "mb-0", "Potential Off-Targets"),
          tags$div(
            class = "d-flex align-items-center gap-2",
            tags$span(class = "small text-muted", "Top 3"),
            tags$div(
              class = "form-check form-switch mb-0",
              tags$input(
                type = "checkbox",
                class = "form-check-input",
                id = ns("show_all_offtargets"),
                onclick = sprintf("Shiny.setInputValue('%s', this.checked, {priority: 'event'})", ns("show_all_offtargets"))
              ),
              tags$label(class = "form-check-label small", `for` = ns("show_all_offtargets"), "Show all")
            )
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

    # Download handler - exports selected probes (or all if none selected)
    output$download_probes <- downloadHandler(
      filename = function() {
        goal <- values$wizard_goal %||% "probes"
        target <- switch(goal,
          "amino_acid" = values$wizard_selection$amino_acid,
          "isoacceptor" = values$wizard_selection$anticodon,
          "specific" = "selected",
          "probes"
        )
        paste0("truenorth_", target, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        probes <- values$wizard_probes
        selected <- values$selected_probes %||% integer()
        coverage_data <- values$coverage_data

        # Filter to selected probes if any are selected
        if (length(selected) > 0) {
          probes <- probes[probes$rank %in% selected, ]
        }

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
          modification_penalty = if ("modification_penalty" %in% names(probes)) round(probes$modification_penalty, 1) else NA,
          overlaps_anticodon = probes$overlaps_anticodon,
          reference_id = probes$reference_id,
          stringsAsFactors = FALSE
        )

        # Add coverage columns if available
        if (!is.null(coverage_data) && "n_targets_hit" %in% names(probes)) {
          export_df$targets_covered_count <- probes$n_targets_hit
          export_df$coverage_percent <- round(probes$coverage_pct, 1)
          export_df$targets_covered <- probes$targets_covered
        }

        # Write with coverage summary header
        con <- file(file, "w")
        on.exit(close(con))

        # Write summary header as comments
        if (length(selected) > 0 && !is.null(coverage_data)) {
          coverage <- get_cumulative_coverage(selected, coverage_data$matrix, coverage_data$target_ids)
          writeLines(sprintf("# TRUENORTH Probe Export - %s", Sys.Date()), con)
          writeLines(sprintf("# Selected probes: %d", length(selected)), con)
          writeLines(sprintf("# Target coverage: %d/%d (%.0f%%)",
                             coverage$n_covered, coverage_data$n_targets, coverage$coverage_pct), con)
          writeLines(sprintf("# Covered targets: %s", paste(coverage$covered_ids, collapse = ", ")), con)
          if (coverage$n_uncovered > 0) {
            writeLines(sprintf("# Uncovered targets: %s", paste(coverage$uncovered_ids, collapse = ", ")), con)
          }
          writeLines("#", con)
        }

        write.csv(export_df, con, row.names = FALSE)
      }
    )
  })
}
