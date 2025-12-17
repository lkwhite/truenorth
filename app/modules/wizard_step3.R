# wizard_step3.R
# Step 3: Feasibility Analysis - show analysis before designing probes

wizardStep3UI <- function(id) {
  ns <- NS(id)

  tagList(
    # Selection context
    uiOutput(ns("selection_context")),

    tags$div(
      class = "text-center mb-4",
      tags$h3("Feasibility Analysis"),
      tags$p(class = "text-muted",
             "Review the analysis before designing probes")
    ),

    uiOutput(ns("feasibility_content"))
  )
}

wizardStep3Server <- function(id, values) {
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

    output$feasibility_content <- renderUI({
      feas <- values$wizard_feasibility

      if (is.null(feas)) {
        return(tags$div(
          class = "text-center text-muted",
          tags$p("Analyzing..."),
          tags$div(class = "spinner-border", role = "status")
        ))
      }

      # Status badge
      status_class <- paste0("feasibility-", feas$status$status)

      tagList(
        # Overall status card
        tags$div(
          class = "card mb-4",
          tags$div(
            class = "card-body text-center",
            tags$h2(
              class = status_class,
              style = "font-weight: bold;",
              feas$status$label
            ),
            tags$p(class = "text-muted mb-0", feas$status$message)
          )
        ),

        # Stats row
        fluidRow(
          # Targets
          column(
            width = 4,
            tags$div(
              class = "card h-100",
              tags$div(
                class = "card-header bg-light",
                tags$strong("Targets")
              ),
              tags$div(
                class = "card-body text-center",
                tags$h2(style = "color: #009E73;", feas$n_targets),
                tags$small(class = "text-muted", "tRNAs to detect")
              )
            )
          ),

          # Non-targets
          column(
            width = 4,
            tags$div(
              class = "card h-100",
              tags$div(
                class = "card-header bg-light",
                tags$strong("Non-targets")
              ),
              tags$div(
                class = "card-body text-center",
                tags$h2(style = "color: #666;", feas$n_nontargets),
                tags$small(class = "text-muted", "tRNAs to avoid")
              )
            )
          ),

          # Conservation
          column(
            width = 4,
            tags$div(
              class = "card h-100",
              tags$div(
                class = "card-header bg-light",
                tags$strong("Target Conservation")
              ),
              tags$div(
                class = "card-body text-center",
                tags$h2(
                  style = paste0("color: ", get_conservation_color(feas$conservation$mean_identity / 100), ";"),
                  sprintf("%.0f%%", feas$conservation$mean_identity %||% 0)
                ),
                tags$small(class = "text-muted", "mean identity")
              )
            )
          )
        ),

        tags$hr(),

        # Detailed analysis
        fluidRow(
          # Specificity info
          column(
            width = 6,
            tags$div(
              class = "card h-100",
              tags$div(
                class = "card-header bg-light",
                icon("shield-alt"), " Specificity"
              ),
              tags$div(
                class = "card-body",
                if (is.null(feas$divergence$n_avoid) || feas$divergence$n_avoid > 0) {
                  # Get the max similarity to closest non-target (lower = better)
                  max_sim <- feas$divergence$max_cross_identity %||% 0
                  tagList(
                    tags$p(
                      tags$strong("Closest non-target similarity: "),
                      tags$span(
                        style = paste0("color: ", get_divergence_color_similarity(max_sim), ";"),
                        sprintf("%.1f%%", max_sim)
                      )
                    ),
                    tags$p(
                      tags$strong("Mean similarity to non-targets: "),
                      sprintf("%.1f%%", feas$divergence$mean_cross_identity %||% 0)
                    ),
                    if (!is.null(feas$divergence$specificity_challenge)) {
                      tags$p(
                        class = "text-muted small",
                        feas$divergence$specificity_challenge
                      )
                    }
                  )
                } else {
                  tags$p(class = "text-success",
                         icon("check-circle"),
                         " No cross-reactivity concerns - all tRNAs are targets")
                }
              )
            )
          ),

          # Best regions
          column(
            width = 6,
            tags$div(
              class = "card h-100",
              tags$div(
                class = "card-header bg-light",
                icon("map-marker-alt"), " Best Probe Regions"
              ),
              tags$div(
                class = "card-body",
                if (!is.null(feas$regions) && nrow(feas$regions) > 0) {
                  top_regions <- head(feas$regions, 3)
                  tags$ul(
                    class = "list-unstyled mb-0",
                    lapply(seq_len(nrow(top_regions)), function(i) {
                      r <- top_regions[i, ]
                      tags$li(
                        class = "mb-2",
                        tags$span(
                          class = "badge bg-primary me-2",
                          paste0("Position ", r$start, "-", r$end)
                        ),
                        tags$small(
                          sprintf("Score: %.2f", r$score)
                        )
                      )
                    })
                  )
                } else {
                  tags$p(class = "text-muted",
                         "Region analysis will be performed during probe design")
                }
              )
            )
          )
        ),

        # Warning messages if needed
        if (feas$status$status %in% c("challenging", "not-feasible")) {
          tags$div(
            class = "alert alert-warning mt-4",
            icon("exclamation-triangle"),
            tags$strong(" Note: "),
            switch(feas$status$status,
              "challenging" = "Probe design may be difficult. Consider selecting fewer targets or a different goal.",
              "not-feasible" = "The selected targets are too similar to non-targets. Consider adjusting your selection."
            )
          )
        },

        # Coverage estimate card
        tags$div(
          class = "card mt-4 border-primary",
          tags$div(
            class = "card-header bg-primary text-white",
            icon("vials"), " Coverage Estimate"
          ),
          tags$div(
            class = "card-body",
            {
              conservation <- feas$conservation$mean_identity %||% 0
              n_targets <- feas$n_targets

              # Estimate probes needed based on conservation
              if (n_targets == 1) {
                estimate_text <- "1 probe"
                estimate_detail <- "Single target - one probe will cover it"
                estimate_class <- "text-success"
              } else if (conservation >= 95) {
                estimate_text <- "1 probe"
                estimate_detail <- sprintf("High conservation (%.0f%%) - one probe should cover all %d targets", conservation, n_targets)
                estimate_class <- "text-success"
              } else if (conservation >= 85) {
                estimate_text <- "1-2 probes"
                estimate_detail <- sprintf("Good conservation (%.0f%%) - may need 1-2 probes for full coverage", conservation)
                estimate_class <- "text-info"
              } else if (conservation >= 70) {
                estimate_text <- "2-3 probes"
                estimate_detail <- sprintf("Moderate conservation (%.0f%%) - expect to need 2-3 probes", conservation)
                estimate_class <- "text-warning"
              } else {
                estimate_text <- "3+ probes"
                estimate_detail <- sprintf("Low conservation (%.0f%%) - multiple probes likely needed", conservation)
                estimate_class <- "text-danger"
              }

              tagList(
                tags$div(
                  class = "text-center mb-3",
                  tags$span(class = paste("fs-4 fw-bold", estimate_class), estimate_text),
                  tags$span(class = "text-muted", " estimated for full coverage")
                ),
                tags$p(class = "text-muted small mb-0", estimate_detail),
                tags$hr(),
                tags$p(class = "small mb-0",
                       icon("info-circle"), " ",
                       "Exact coverage will be calculated after probe design. ",
                       "In the next step, you can select which probes to include and see actual target coverage."
                )
              )
            }
          )
        ),

        # Design parameters
        tags$div(
          class = "card mt-4",
          tags$div(
            class = "card-header bg-light",
            icon("cog"), " Design Parameters"
          ),
          tags$div(
            class = "card-body",
            fluidRow(
              column(
                width = 4,
                numericInput(ns("probe_length"), "Probe length (nt)", value = 20, min = 15, max = 30)
              ),
              column(
                width = 4,
                selectInput(ns("region_pref"), "Target region",
                            choices = c("Any region" = "any",
                                        "5' half" = "5prime",
                                        "3' half" = "3prime"),
                            selected = "any")
              ),
              column(
                width = 4,
                tags$label(class = "form-label", "Anticodon avoidance"),
                tags$div(
                  class = "form-check form-switch",
                  tags$input(type = "checkbox", class = "form-check-input",
                             id = ns("avoid_anticodon"), checked = "checked"),
                  tags$label(class = "form-check-label", `for` = ns("avoid_anticodon"),
                             "Avoid anticodon region")
                ),
                tags$small(class = "text-muted d-block mt-1",
                           "Heavily modified, poor hybridization")
              )
            )
          )
        )
      )
    })

    # Expose design parameters
    observe({
      # Store parameters for design step
      values$design_params <- list(
        probe_length = input$probe_length %||% 20,
        region_pref = input$region_pref %||% "any",
        avoid_anticodon = if (!is.null(input$avoid_anticodon)) input$avoid_anticodon else TRUE
      )
    })
  })
}

# Helper functions for color coding
get_conservation_color <- function(conservation) {
  if (is.null(conservation) || is.na(conservation)) return("#666")
  if (conservation >= 0.95) return("#009E73")  # Excellent
  if (conservation >= 0.85) return("#56B4E9")  # Good
  if (conservation >= 0.70) return("#E69F00")  # Moderate
  return("#D55E00")  # Low
}

# For divergence, higher similarity to non-targets is BAD (less specific)
get_divergence_color_similarity <- function(similarity) {
  if (is.null(similarity) || is.na(similarity)) return("#666")
  # Lower similarity to non-targets = better specificity
  if (similarity <= 70) return("#009E73")   # Excellent specificity
  if (similarity <= 80) return("#56B4E9")   # Good specificity
  if (similarity <= 85) return("#E69F00")   # Moderate specificity
  if (similarity <= 90) return("#D55E00")   # Challenging
  return("#CC79A7")  # Not feasible - too similar
}
