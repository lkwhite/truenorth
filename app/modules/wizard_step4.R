# wizard_step4.R
# Step 4: Design Results - probe table with coverage selection and export

wizardStep4UI <- function(id) {
  ns <- NS(id)

  tagList(
    # Help section
    collapsible_help(
      ns("help_step4"),
      "Tips for selecting probes",
      help_step4_ui(),
      collapsed = TRUE
    ),

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

    # Track probe display mode: "all" or "best"
    probe_display_mode <- reactiveVal("best")

    # Handle All/Best toggle buttons
    observeEvent(input$show_all, {
      probe_display_mode("all")
    })

    observeEvent(input$show_best, {
      probe_display_mode("best")
    })

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
      distinguish_mode <- values$distinguish_mode %||% FALSE

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

      # DISTINGUISH MODE: Show probes grouped by target
      if (distinguish_mode && !is.null(values$probe_sets)) {
        return(render_distinguish_results(ns, probes, values$probe_sets, values$trna_data))
      }

      # TOGETHER MODE: Standard display with coverage tracking
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
            uiOutput(ns("coverage_tracker")),

            # Collapsible coverage matrix
            tags$hr(class = "my-3"),
            tags$div(
              tags$div(
                class = "d-flex justify-content-between align-items-center",
                tags$a(
                  class = "text-decoration-none",
                  `data-bs-toggle` = "collapse",
                  href = paste0("#", ns("coverage_matrix_collapse")),
                  role = "button",
                  `aria-expanded` = "true",
                  tags$span(icon("th"), " Coverage Matrix"),
                  tags$small(class = "text-muted ms-2", "(click to collapse)")
                ),
                tags$div(
                  class = "d-flex gap-2",
                  actionButton(
                    ns("suggest_minimum"),
                    "Suggest minimum set",
                    class = "btn btn-sm btn-outline-success",
                    icon = icon("magic")
                  ),
                  actionButton(
                    ns("clear_selection"),
                    "Clear",
                    class = "btn btn-sm btn-outline-secondary"
                  )
                )
              ),
              tags$div(
                id = ns("coverage_matrix_collapse"),
                class = "collapse show mt-3",
                tags$p(class = "small text-muted mb-2",
                  "Click on probe columns (#1, #2, etc.) to select/deselect probes"
                ),
                uiOutput(ns("coverage_matrix_view"))
              )
            )
          )
        ),

        # Probe table with filter controls
        tags$div(
          class = "card",
          tags$div(
            class = "card-header bg-light",
            tags$div(
              class = "d-flex justify-content-between align-items-center",
              tags$strong("Probe Details"),
              uiOutput(ns("probe_toggle_buttons"))
            ),
            # Filter controls row
            tags$div(
              class = "mt-2 pt-2 border-top d-flex flex-wrap gap-3 align-items-center small",
              # Specificity filter
              tags$div(
                class = "d-flex align-items-center gap-2",
                tags$label(class = "mb-0 text-muted", "Min diff-AA mismatches:"),
                tags$select(
                  id = ns("min_diff_aa_filter"),
                  class = "form-select form-select-sm",
                  style = "width: 70px;",
                  onchange = sprintf("Shiny.setInputValue('%s', this.value, {priority: 'event'})", ns("min_diff_aa_filter")),
                  tags$option(value = "0", "Any"),
                  tags$option(value = "3", "3+"),
                  tags$option(value = "4", "4+"),
                  tags$option(value = "5", "5+", selected = "selected"),
                  tags$option(value = "6", "6+")
                )
              ),
              # Clustering toggle
              tags$div(
                class = "form-check",
                tags$input(
                  type = "checkbox",
                  class = "form-check-input",
                  id = ns("cluster_probes"),
                  onclick = sprintf("Shiny.setInputValue('%s', this.checked, {priority: 'event'})", ns("cluster_probes"))
                ),
                tags$label(
                  class = "form-check-label text-muted",
                  `for` = ns("cluster_probes"),
                  "Deduplicate (unique regions)"
                )
              )
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

    # Toggle buttons for All/Best view
    output$probe_toggle_buttons <- renderUI({
      mode <- probe_display_mode()
      tags$div(
        class = "btn-group btn-group-sm",
        actionButton(
          ns("show_all"), "All",
          class = if (mode == "all") "btn-secondary" else "btn-outline-secondary"
        ),
        actionButton(
          ns("show_best"), "Best Only",
          class = if (mode == "best") "btn-secondary" else "btn-outline-secondary"
        )
      )
    })

    # Probe table with multi-select (click to select for probe set AND view details)
    output$probe_table <- DT::renderDataTable({
      req(values$wizard_probes)
      probes <- values$wizard_probes
      coverage_data <- values$coverage_data

      # Apply specificity filter
      min_diff_aa <- as.integer(input$min_diff_aa_filter %||% "0")
      if (min_diff_aa > 0 && "min_diff_aa_mismatches" %in% names(probes)) {
        probes <- probes[
          is.na(probes$min_diff_aa_mismatches) |
          probes$min_diff_aa_mismatches >= min_diff_aa,
        ]
      }

      # Apply clustering/deduplication
      if (isTRUE(input$cluster_probes) && nrow(probes) > 1) {
        probes <- cluster_probes_by_position(probes, min_overlap = 0.7)
      }

      # Filter based on display mode
      mode <- probe_display_mode()
      if (mode == "best") {
        # Show top 10 probes by rank
        probes <- head(probes[order(probes$rank), ], 10)
      }

      # Handle empty result after filtering
      if (nrow(probes) == 0) {
        return(DT::datatable(
          data.frame(Message = "No probes match the current filters. Try relaxing the specificity threshold."),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }

      # Format coverage columns
      if (!is.null(coverage_data) && "n_targets_hit" %in% names(probes)) {
        n_targets <- coverage_data$n_targets
        covers_col <- sprintf("%d/%d (%.0f%%)", probes$n_targets_hit, n_targets, probes$coverage_pct)
      } else {
        covers_col <- rep("-", nrow(probes))
      }

      # Format Diff-AA off-target column
      if ("min_diff_aa_mismatches" %in% names(probes)) {
        diff_aa_col <- sapply(seq_len(nrow(probes)), function(i) {
          mm <- probes$min_diff_aa_mismatches[i]
          aa <- probes$closest_diff_aa_amino[i]
          if (is.na(mm)) {
            "-"
          } else {
            # Color code: green >=5, yellow 3-4, red <3
            color <- if (mm >= 5) "#009E73" else if (mm >= 3) "#E69F00" else "#D55E00"
            sprintf('<span style="color:%s; font-weight:bold;">%d</span> <small>(%s)</small>',
                    color, mm, if (is.na(aa)) "?" else aa)
          }
        })
      } else {
        diff_aa_col <- rep("-", nrow(probes))
      }

      # Format for display
      display_df <- data.frame(
        Rank = probes$rank,
        Sequence = probes$probe_sequence,
        Region = probes$trna_region,
        Covers = covers_col,
        `Diff-AA` = diff_aa_col,
        Tm = sprintf("%.1f", probes$tm_nn),
        Quality = probes$quality,
        stringsAsFactors = FALSE,
        check.names = FALSE
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
        selection = list(mode = "multiple"),  # Selection managed via proxy
        options = list(
          pageLength = 25,  # Show all typical probe sets on one page
          dom = "ti",  # table + info only, no pagination for small sets
          columnDefs = list(
            list(className = "dt-center", targets = c(0, 3, 4, 5, 6))  # Rank, Covers, Diff-AA, Tm, Quality
          )
        ),
        rownames = FALSE,
        escape = FALSE  # Allow HTML/unicode in cells
      )
    })

    # Create table proxy for updating selection without re-rendering
    probe_table_proxy <- DT::dataTableProxy("probe_table")

    # Sync FROM values$selected_probes TO table (one-way: values -> table)
    observeEvent(values$selected_probes, {
      req(values$wizard_probes)
      selected_ranks <- values$selected_probes %||% integer()
      # Find which rows have these ranks
      selected_rows <- which(values$wizard_probes$rank %in% selected_ranks)
      DT::selectRows(probe_table_proxy, selected_rows)
    }, ignoreInit = TRUE)

    # Sync FROM table TO values (only on explicit table clicks)
    # Use a flag to prevent circular updates
    observeEvent(input$probe_table_rows_selected, {
      req(values$wizard_probes)
      # Only process if this is a user click, not a programmatic update
      # We detect this by checking if the selection actually differs from what we set
      selected_rows <- input$probe_table_rows_selected
      if (is.null(selected_rows)) selected_rows <- integer()

      new_ranks <- if (length(selected_rows) > 0) {
        values$wizard_probes$rank[selected_rows]
      } else {
        integer()
      }

      current_ranks <- values$selected_probes %||% integer()

      # Only update if different (avoids circular updates)
      if (!setequal(new_ranks, current_ranks)) {
        values$selected_probes <- new_ranks
      }
    }, ignoreInit = TRUE)

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

    # Suggest minimum probe set using greedy set cover
    observeEvent(input$suggest_minimum, {
      coverage_data <- values$coverage_data
      probes <- values$wizard_probes

      if (is.null(coverage_data) || is.null(probes)) {
        showNotification("Coverage data not available", type = "warning")
        return()
      }

      # Use individual target matrix and build group coverage fresh
      # (same logic as visualization to ensure consistency)
      ind_matrix <- coverage_data$matrix
      target_sequences <- coverage_data$target_sequences
      target_ids <- coverage_data$target_ids

      if (is.null(ind_matrix) || nrow(ind_matrix) == 0) {
        showNotification("No coverage data to analyze", type = "warning")
        return()
      }

      # Build group coverage matrix using same logic as visualization
      if (!is.null(target_sequences) && length(target_sequences) == length(target_ids)) {
        unique_seqs <- unique(target_sequences)
        n_groups <- length(unique_seqs)
        n_probes <- nrow(ind_matrix)

        # Build group membership (which columns belong to each group)
        group_cols <- lapply(unique_seqs, function(s) {
          which(target_sequences == s)
        })

        # Build group coverage matrix
        coverage_matrix <- matrix(FALSE, nrow = n_probes, ncol = n_groups)
        for (g in seq_len(n_groups)) {
          cols <- group_cols[[g]]
          for (p in seq_len(n_probes)) {
            coverage_matrix[p, g] <- any(ind_matrix[p, cols])
          }
        }
      } else {
        # No grouping - use individual matrix
        coverage_matrix <- ind_matrix
        n_groups <- ncol(coverage_matrix)
        unique_seqs <- target_ids
      }

      # Greedy set cover algorithm
      n_probes <- nrow(coverage_matrix)
      n_targets <- ncol(coverage_matrix)
      selected_indices <- integer()
      covered <- rep(FALSE, n_targets)

      while (!all(covered) && length(selected_indices) < n_probes) {
        # Find probe that covers the most uncovered targets
        uncovered_coverage <- sapply(seq_len(n_probes), function(p) {
          if (p %in% selected_indices) return(-1)  # Already selected
          sum(coverage_matrix[p, ] & !covered)
        })

        best_probe <- which.max(uncovered_coverage)
        best_new <- uncovered_coverage[best_probe]

        if (best_new == 0) break  # No probe can cover any more targets

        selected_indices <- c(selected_indices, best_probe)
        covered <- covered | coverage_matrix[best_probe, ]
      }

      # Convert row indices to probe ranks
      selected_ranks <- probes$rank[selected_indices]

      # Update selection
      values$selected_probes <- sort(selected_ranks)

      # Notify user
      n_covered <- sum(covered)
      n_total <- n_targets
      if (all(covered)) {
        showNotification(
          sprintf("Selected %d probes for full coverage of all %d target groups",
                  length(selected_ranks), n_total),
          type = "message"
        )
      } else {
        showNotification(
          sprintf("Selected %d probes covering %d/%d target groups (some targets may be uncoverable)",
                  length(selected_ranks), n_covered, n_total),
          type = "warning"
        )
      }
    })

    # Coverage matrix visualization
    output$coverage_matrix_view <- renderUI({
      coverage_data <- values$coverage_data
      probes <- values$wizard_probes
      selected <- values$selected_probes %||% integer()

      if (is.null(coverage_data) || is.null(probes)) {
        return(tags$p(class = "text-muted", "Coverage data not available"))
      }

      render_coverage_matrix(
        coverage_matrix = coverage_data$matrix,
        probe_ranks = probes$rank,
        target_ids = coverage_data$target_ids,
        target_sequences = coverage_data$target_sequences,
        selected_probes = selected,
        max_probes = 15,
        max_targets = 50,  # Show all groups (sequences are collapsed)
        ns = ns  # Pass namespace for clickable probe headers
      )
    })

    # Handle clicks on probe columns in coverage matrix
    observeEvent(input$matrix_probe_click, {
      click_data <- input$matrix_probe_click
      if (is.null(click_data)) return()

      rank <- click_data$rank
      current_selected <- values$selected_probes %||% integer()

      # Toggle selection
      if (rank %in% current_selected) {
        values$selected_probes <- setdiff(current_selected, rank)
      } else {
        values$selected_probes <- c(current_selected, rank)
      }
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

      # Use last clicked or first selected for detail view
      detail_row <- input$probe_table_row_last_clicked
      if (is.null(detail_row) || !detail_row %in% selected_row) {
        detail_row <- selected_row[1]
      }
      probe <- values$wizard_probes[detail_row, ]
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

    # Reactive off-target diagram with toggle - now separated by amino acid
    output$offtarget_diagram <- renderUI({
      req(values$wizard_probes)

      selected_row <- input$probe_table_rows_selected
      if (is.null(selected_row) || length(selected_row) == 0) {
        return(NULL)
      }

      # Use last clicked or first selected for detail view
      detail_row <- input$probe_table_row_last_clicked
      if (is.null(detail_row) || !detail_row %in% selected_row) {
        detail_row <- selected_row[1]
      }
      probe <- values$wizard_probes[detail_row, ]
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

      # Get target amino acids
      target_trnas <- all_trnas[all_trnas$id %in% target_ids, ]
      target_amino_acids <- unique(target_trnas$amino_acid)

      # Separate same-AA and different-AA off-targets
      same_aa_df <- off_target_df[off_target_df$amino_acid %in% target_amino_acids, ]
      diff_aa_df <- off_target_df[!off_target_df$amino_acid %in% target_amino_acids, ]

      # Calculate mismatches for each off-target
      probe_reversed <- paste(rev(strsplit(probe$probe_sequence, "")[[1]]), collapse = "")
      probe_chars <- strsplit(probe_reversed, "")[[1]]

      calc_mismatches <- function(df) {
        if (nrow(df) == 0) return(df)
        df$binding_region <- substr(df$sequence, probe$start, probe$end)
        df$n_mismatches <- sapply(df$binding_region, function(region) {
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
        df[order(df$n_mismatches), ]
      }

      same_aa_df <- calc_mismatches(same_aa_df)
      diff_aa_df <- calc_mismatches(diff_aa_df)

      # Use toggle to determine max_show
      show_all <- input$show_all_offtargets %||% FALSE
      max_show <- if (show_all) Inf else 3

      # Build output sections
      sections <- list()

      # CRITICAL: Different amino acid off-targets (must maintain specificity)
      if (nrow(diff_aa_df) > 0) {
        closest_diff <- if (show_all) diff_aa_df else head(diff_aa_df, max_show)
        min_diff_mm <- min(closest_diff$n_mismatches)

        diff_alert <- if (min_diff_mm >= 5) {
          tags$div(
            class = "alert alert-success py-2 mb-2",
            icon("check-circle"),
            sprintf(" Good specificity: %d+ mismatches to other amino acids", min_diff_mm)
          )
        } else if (min_diff_mm >= 3) {
          tags$div(
            class = "alert alert-warning py-2 mb-2",
            icon("exclamation-triangle"),
            sprintf(" Moderate: %d mismatches to closest different amino acid", min_diff_mm)
          )
        } else {
          tags$div(
            class = "alert alert-danger py-2 mb-2",
            icon("exclamation-circle"),
            sprintf(" Warning: Only %d mismatch(es) to a different amino acid - may cross-react!", min_diff_mm)
          )
        }

        sections$diff_aa <- tagList(
          tags$h6(
            class = "mt-3 mb-2 text-danger",
            icon("exclamation-triangle", class = "me-1"),
            sprintf("Different Amino Acid Off-targets (%d total)", nrow(diff_aa_df)),
            tags$small(class = "text-muted ms-2", "- must maintain specificity")
          ),
          diff_alert,
          render_multi_target_hybridization(
            targets_df = closest_diff,
            probe_sequence = probe$probe_sequence,
            start = probe$start,
            end = probe$end,
            reference_id = NULL,
            max_show = max_show
          )
        )
      }

      # ALLOWABLE: Same amino acid off-targets (covered by other probes in set)
      if (nrow(same_aa_df) > 0) {
        closest_same <- if (show_all) same_aa_df else head(same_aa_df, max_show)
        min_same_mm <- min(closest_same$n_mismatches)

        same_info <- if (min_same_mm <= 2) {
          tags$div(
            class = "alert alert-info py-2 mb-2",
            icon("info-circle"),
            sprintf(" Closest same-AA off-target has %d mismatch(es) - may need additional probe for coverage", min_same_mm)
          )
        } else {
          tags$div(
            class = "alert alert-light py-2 mb-2 border",
            icon("check"),
            sprintf(" Same-AA off-targets have %d+ mismatches - use other probes in set to cover", min_same_mm)
          )
        }

        sections$same_aa <- tagList(
          tags$h6(
            class = "mt-4 mb-2 text-info",
            icon("layer-group", class = "me-1"),
            sprintf("Same Amino Acid Off-targets (%d total)", nrow(same_aa_df)),
            tags$small(class = "text-muted ms-2", "- allowable, cover with other probes")
          ),
          same_info,
          render_multi_target_hybridization(
            targets_df = closest_same,
            probe_sequence = probe$probe_sequence,
            start = probe$start,
            end = probe$end,
            reference_id = NULL,
            max_show = max_show
          )
        )
      }

      # Return both sections (different-AA first as it's more critical)
      if (length(sections) == 0) {
        return(tags$div(
          class = "alert alert-success py-2",
          icon("check-circle"),
          " No off-targets found"
        ))
      }

      tagList(sections$diff_aa, sections$same_aa)
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

# =============================================================================
# Helper function for distinguish mode results
# =============================================================================

render_distinguish_results <- function(ns, all_probes, probe_sets, trna_data) {
  n_targets <- length(probe_sets)
  n_total_probes <- nrow(all_probes)
  target_ids <- names(probe_sets)

  tagList(
    # CSS for clickable probe rows
    tags$style(HTML("
      .distinguish-probe-row {
        cursor: pointer;
        transition: background-color 0.15s;
      }
      .distinguish-probe-row:hover {
        background-color: #e3f2fd !important;
      }
      .distinguish-probe-row.selected {
        background-color: #bbdefb !important;
      }
      .distinguish-detail {
        display: none;
        animation: fadeIn 0.2s;
      }
      .distinguish-detail.show {
        display: block;
      }
      @keyframes fadeIn {
        from { opacity: 0; }
        to { opacity: 1; }
      }
    ")),

    # JavaScript for row selection
    tags$script(HTML("
      function selectDistinguishProbe(groupId, rowIdx) {
        // Deselect all rows in this group
        document.querySelectorAll('#group-' + groupId + ' .distinguish-probe-row').forEach(function(row) {
          row.classList.remove('selected');
        });
        // Select clicked row
        document.getElementById('probe-row-' + groupId + '-' + rowIdx).classList.add('selected');
        // Hide all details in this group
        document.querySelectorAll('#group-' + groupId + ' .distinguish-detail').forEach(function(detail) {
          detail.classList.remove('show');
        });
        // Show selected detail
        document.getElementById('probe-detail-' + groupId + '-' + rowIdx).classList.add('show');
      }
    ")),

    # Summary header
    tags$div(
      class = "alert alert-info mb-4",
      tags$div(
        class = "d-flex align-items-center",
        icon("code-branch", class = "me-2"),
        tags$div(
          tags$strong("Distinguish Mode"),
          tags$span(class = "ms-2",
            sprintf("Designed %d probe sets to distinguish between %d isodecoders",
                    n_targets, n_targets))
        )
      )
    ),

    # Export button
    tags$div(
      class = "mb-4 text-end",
      downloadButton(ns("download_probes"), "Export All Probes (CSV)",
                     class = "btn-primary btn-sm")
    ),

    # Probe sets grouped by target
    lapply(seq_along(target_ids), function(idx) {
      target_id <- target_ids[idx]
      probes <- probe_sets[[target_id]]
      if (is.null(probes) || nrow(probes) == 0) return(NULL)

      target_label <- format_trna_id(target_id)
      n_probes <- nrow(probes)
      other_targets <- setdiff(target_ids, target_id)
      other_targets_df <- trna_data[trna_data$id %in% other_targets, ]

      tags$div(
        id = paste0("group-", idx),
        class = "card mb-4",
        tags$div(
          class = "card-header bg-primary text-white",
          tags$div(
            class = "d-flex justify-content-between align-items-center",
            tags$span(
              icon("crosshairs", class = "me-2"),
              tags$strong(target_label),
              tags$small(class = "ms-2 opacity-75",
                sprintf("(distinguishes from %s)", paste(sapply(other_targets, format_trna_id), collapse = ", ")))
            ),
            tags$span(class = "badge bg-light text-primary",
                      paste(n_probes, "probes"))
          )
        ),
        tags$div(
          class = "card-body",
          # Probe table as HTML
          tags$table(
            class = "table table-sm table-hover",
            tags$thead(
              tags$tr(
                tags$th("Rank", style = "width: 60px;"),
                tags$th("Sequence"),
                tags$th("Region", style = "width: 100px;"),
                tags$th("Tm", style = "width: 70px;"),
                tags$th("GC", style = "width: 60px;"),
                tags$th("Quality", style = "width: 80px;")
              )
            ),
            tags$tbody(
              lapply(seq_len(nrow(probes)), function(row_idx) {
                probe <- probes[row_idx, ]
                region_display <- probe$trna_region
                if ("modification_penalty" %in% names(probe)) {
                  if (probe$modification_penalty >= 20) {
                    region_display <- paste0(region_display, " \u26A0\u26A0")
                  } else if (probe$modification_penalty >= 10) {
                    region_display <- paste0(region_display, " \u26A0")
                  }
                }
                tags$tr(
                  id = paste0("probe-row-", idx, "-", row_idx),
                  class = "distinguish-probe-row",
                  onclick = sprintf("selectDistinguishProbe(%d, %d)", idx, row_idx),
                  tags$td(probe$rank),
                  tags$td(tags$code(probe$probe_sequence, style = "font-size: 0.9em;")),
                  tags$td(region_display),
                  tags$td(sprintf("%.1f", probe$tm_nn)),
                  tags$td(sprintf("%.0f%%", probe$gc_content)),
                  tags$td(probe$quality)
                )
              })
            )
          ),

          # Instruction text
          tags$p(class = "text-muted small mt-2 mb-3",
                 icon("hand-pointer"), " Click a probe row to see cross-hybridization details"),

          # Detail views (hidden by default, shown when row clicked)
          lapply(seq_len(nrow(probes)), function(row_idx) {
            probe <- probes[row_idx, ]
            tags$div(
              id = paste0("probe-detail-", idx, "-", row_idx),
              class = "distinguish-detail border rounded p-3 bg-light",
              render_distinguish_cross_hyb(probe, target_id, other_targets_df)
            )
          })
        )
      )
    })
  )
}

# =============================================================================
# Helper function for distinguish mode cross-hybridization display
# =============================================================================

#' Render cross-hybridization analysis for distinguish mode
#'
#' Shows how a probe designed for one target would bind to other selected targets
#'
#' @param probe Single probe row (data frame)
#' @param target_id The target this probe was designed for
#' @param other_targets_df Data frame of other selected targets
#' @return HTML tags showing cross-hybridization potential
render_distinguish_cross_hyb <- function(probe, target_id, other_targets_df) {
  if (is.null(other_targets_df) || nrow(other_targets_df) == 0) {
    return(tags$p(class = "text-muted", "No other targets to analyze"))
  }

  # Probe binds antiparallel - reverse for comparison
  probe_reversed <- paste(rev(strsplit(probe$probe_sequence, "")[[1]]), collapse = "")
  probe_chars <- strsplit(probe_reversed, "")[[1]]
  probe_len <- length(probe_chars)

  # Calculate base Tm for this probe
  base_tm <- probe$tm_nn

  # Extract binding regions and calculate mismatches
  other_targets_df$binding_region <- substr(other_targets_df$sequence, probe$start, probe$end)

  other_targets_df$n_mismatches <- sapply(other_targets_df$binding_region, function(region) {
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

  # Estimate Tm with mismatches (~5°C reduction per mismatch)
  other_targets_df$estimated_tm <- base_tm - (other_targets_df$n_mismatches * 5)

  # Sort by mismatches
  other_targets_df <- other_targets_df[order(other_targets_df$n_mismatches), ]

  # Build visualization for each other target
  target_elements <- lapply(seq_len(nrow(other_targets_df)), function(i) {
    target <- other_targets_df[i, ]
    n_mm <- target$n_mismatches
    est_tm <- target$estimated_tm
    region <- target$binding_region
    region_chars <- strsplit(region, "")[[1]]

    # Determine specificity level
    specificity <- if (n_mm >= 5) {
      list(class = "text-success", icon = "check-circle",
           label = "Good specificity", desc = "Unlikely to cross-hybridize")
    } else if (n_mm >= 3) {
      list(class = "text-warning", icon = "exclamation-triangle",
           label = "Moderate", desc = "May show some cross-hybridization")
    } else {
      list(class = "text-danger", icon = "exclamation-circle",
           label = "Poor specificity", desc = "Will likely cross-hybridize")
    }

    # Build sequence with mismatches highlighted and track positions
    mismatch_positions <- c()
    seq_html <- lapply(seq_along(region_chars), function(j) {
      target_base <- toupper(region_chars[j])
      probe_base <- if (j <= length(probe_chars)) toupper(probe_chars[j]) else "?"

      is_match <- (target_base == "A" && probe_base == "T") ||
                  (target_base == "T" && probe_base == "A") ||
                  (target_base == "G" && probe_base == "C") ||
                  (target_base == "C" && probe_base == "G")

      if (is_match) {
        tags$span(style = "color: #009E73;", target_base)  # Match - green
      } else {
        tags$span(
          style = "color: #D55E00; font-weight: bold; background: #FFF3CD;",
          title = sprintf("Position %d: %s in target vs %s in probe", j, target_base, probe_base),
          target_base
        )  # Mismatch - red with highlight
      }
    })

    # Find mismatch positions for display
    mm_pos <- which(sapply(seq_along(region_chars), function(j) {
      target_base <- toupper(region_chars[j])
      probe_base <- if (j <= length(probe_chars)) toupper(probe_chars[j]) else "?"
      !((target_base == "A" && probe_base == "T") ||
        (target_base == "T" && probe_base == "A") ||
        (target_base == "G" && probe_base == "C") ||
        (target_base == "C" && probe_base == "G"))
    }))

    tags$div(
      class = "mb-3 p-2 border rounded",
      style = "background: #f8f9fa;",
      # Header row
      tags$div(
        class = "d-flex justify-content-between align-items-center mb-2",
        tags$span(
          tags$strong(format_trna_id(target$id)),
          tags$span(
            class = "badge bg-secondary ms-2",
            sprintf("%d mismatch%s", n_mm, if (n_mm == 1) "" else "es")
          )
        ),
        tags$span(
          class = paste("small", specificity$class),
          icon(specificity$icon, class = "me-1"),
          specificity$label
        )
      ),
      # Stats row
      tags$div(
        class = "small text-muted mb-2",
        tags$span(
          class = "me-3",
          icon("thermometer-half", class = "me-1"),
          sprintf("Est. binding Tm: %.0f°C", est_tm),
          if (est_tm < 37) tags$span(class = "text-success ms-1", "(below 37°C - minimal binding)") else NULL
        ),
        if (length(mm_pos) > 0 && length(mm_pos) <= 6) {
          tags$span(
            icon("map-marker-alt", class = "me-1"),
            sprintf("Mismatches at positions: %s", paste(mm_pos, collapse = ", "))
          )
        }
      ),
      # Sequence alignment view
      tags$div(
        style = "font-family: monospace; font-size: 0.85em; background: white; padding: 8px; border-radius: 4px;",
        tags$div(
          tags$span(class = "text-muted", "Target 5' "),
          tags$span(seq_html),
          tags$span(class = "text-muted", " 3'")
        ),
        tags$div(
          class = "text-muted",
          tags$span("Probe  3' "),
          tags$span(probe_reversed),
          tags$span(" 5'")
        )
      )
    )
  })

  # Summary header
  min_mm <- min(other_targets_df$n_mismatches)
  min_est_tm <- min(other_targets_df$estimated_tm)
  summary_class <- if (min_mm >= 5) "alert-success" else if (min_mm >= 3) "alert-warning" else "alert-danger"
  summary_icon <- if (min_mm >= 5) "check-circle" else "exclamation-triangle"

  tagList(
    # Probe properties header
    tags$div(
      class = "mb-3 p-2 border-start border-primary border-3",
      style = "background: #e7f1ff;",
      tags$div(
        class = "d-flex flex-wrap gap-3",
        tags$div(
          tags$small(class = "text-muted d-block", "Probe Sequence"),
          tags$code(probe$probe_sequence, style = "font-size: 1.05em;")
        ),
        tags$div(
          tags$small(class = "text-muted d-block", "Region"),
          tags$strong(probe$trna_region)
        ),
        tags$div(
          tags$small(class = "text-muted d-block", "Position"),
          tags$strong(sprintf("%d-%d", probe$start, probe$end))
        ),
        tags$div(
          tags$small(class = "text-muted d-block", "Tm (perfect match)"),
          tags$strong(sprintf("%.1f°C", probe$tm_nn))
        ),
        tags$div(
          tags$small(class = "text-muted d-block", "GC Content"),
          tags$strong(sprintf("%.0f%%", probe$gc_content))
        )
      )
    ),

    # Summary alert
    tags$div(
      class = paste("alert py-2 mb-3", summary_class),
      icon(summary_icon, class = "me-1"),
      if (min_mm >= 5) {
        sprintf("This probe has %d+ mismatches with all other selected targets - good for distinguishing", min_mm)
      } else if (min_mm >= 3) {
        sprintf("Closest other target has %d mismatches - moderate differentiation", min_mm)
      } else {
        sprintf("Warning: Closest other target has only %d mismatch(es) - may cross-hybridize", min_mm)
      },
      tags$br(),
      tags$small(
        class = "opacity-75",
        sprintf("Estimated off-target Tm range: %.0f°C to %.0f°C",
                min(other_targets_df$estimated_tm), max(other_targets_df$estimated_tm))
      )
    ),

    # Cross-hybridization with other targets
    tags$h6(class = "mb-2",
            icon("dna", class = "me-1"),
            sprintf("Cross-hybridization with %d other selected target%s:",
                    nrow(other_targets_df), if (nrow(other_targets_df) == 1) "" else "s")),
    tags$div(
      style = "max-height: 350px; overflow-y: auto;",
      target_elements
    )
  )
}
