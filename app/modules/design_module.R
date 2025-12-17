# design_module.R
# Shiny module for probe design parameters

# =============================================================================
# UI
# =============================================================================

designUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # Left column: Selection summary and parameters
      column(6,
        # Selection summary
        tags$div(
          class = "card mb-3",
          tags$div(
            class = "card-header",
            "Target Selection"
          ),
          tags$div(
            class = "card-body",
            uiOutput(ns("selection_summary"))
          )
        ),

        # Probe parameters
        tags$div(
          class = "card",
          tags$div(
            class = "card-header",
            "Probe Parameters"
          ),
          tags$div(
            class = "card-body",

            # Probe length
            sliderInput(
              ns("probe_length"),
              "Probe Length (nt)",
              min = 15,
              max = 30,
              value = c(20, 25),
              step = 1
            ),

            # Target region
            selectInput(
              ns("target_region"),
              "Target Region",
              choices = c(
                "Entire sequence" = "all",
                "5' end (positions 1-25)" = "5prime",
                "3' end (last 25 positions)" = "3prime",
                "Middle region" = "middle",
                "Divergent regions (recommended)" = "divergent"
              ),
              selected = "divergent"
            ),

            # Tm range
            sliderInput(
              ns("tm_range"),
              "Tm Range (°C)",
              min = 40,
              max = 75,
              value = c(50, 65),
              step = 1
            ),

            # Number of probes
            numericInput(
              ns("n_probes"),
              "Number of probes to return",
              value = 10,
              min = 1,
              max = 50,
              step = 1
            ),

            # Design button
            tags$div(
              class = "d-grid gap-2 mt-4",
              actionButton(
                ns("design"),
                "Design Probes",
                class = "btn-primary btn-lg",
                icon = icon("flask")
              )
            )
          )
        )
      ),

      # Right column: Help and tips
      column(6,
        tags$div(
          class = "card",
          tags$div(
            class = "card-header",
            "Design Tips"
          ),
          tags$div(
            class = "card-body",
            tags$h6("Probe Length"),
            tags$p(
              "20-25 nt is typical for tRNA northern blots. ",
              "Shorter probes have lower Tm, longer probes may have more off-targets."
            ),

            tags$h6("Target Region"),
            tags$p(
              tags$strong("Divergent regions"), " (recommended): ",
              "Automatically finds regions that differ most between your desired ",
              "targets and those you want to avoid. Best for isoacceptor discrimination."
            ),
            tags$p(
              tags$strong("5' or 3' end"), ": ",
              "Useful if you know the region you want to target. ",
              "The 3' end includes the CCA tail."
            ),

            tags$h6("Tm Range"),
            tags$p(
              "50-65°C is typical for northern blots. ",
              "Higher Tm probes bind more strongly but may have more off-targets."
            ),

            tags$hr(),
            tags$p(
              class = "text-muted small",
              "Tm is calculated using the nearest-neighbor method (SantaLucia 1998) ",
              "assuming 50 nM probe and 50 mM Na+."
            )
          )
        )
      )
    )
  )
}

# =============================================================================
# Server
# =============================================================================

designServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Local reactive for probes
    local <- reactiveValues(
      probes = NULL
    )

    # -------------------------------------------------------------------------
    # Selection summary
    # -------------------------------------------------------------------------
    output$selection_summary <- renderUI({
      n_desired <- sum(unlist(values$selection_types) == "desired")
      n_avoid <- sum(unlist(values$selection_types) == "avoid")

      if (n_desired == 0 && n_avoid == 0) {
        return(tags$div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          " No targets selected. Go to ",
          tags$strong("Browse tRNAs"),
          " to select which tRNAs you want to hit (desired) and avoid."
        ))
      }

      # Get details of selection
      desired_ids <- names(values$selection_types)[unlist(values$selection_types) == "desired"]
      avoid_ids <- names(values$selection_types)[unlist(values$selection_types) == "avoid"]

      tags$div(
        tags$div(
          style = paste0("color: ", COLORS$desired, "; margin-bottom: 10px;"),
          icon("check"),
          tags$strong(n_desired, " desired target(s)"),
          if (n_desired > 0 && n_desired <= 5) {
            tags$ul(
              class = "small mb-0",
              lapply(desired_ids, function(x) tags$li(x))
            )
          } else if (n_desired > 5) {
            tags$span(class = "small", paste0(" (", desired_ids[1], ", ", desired_ids[2], ", ...)"))
          }
        ),
        tags$div(
          style = paste0("color: ", COLORS$avoid, ";"),
          icon("times"),
          tags$strong(n_avoid, " target(s) to avoid"),
          if (n_avoid > 0 && n_avoid <= 5) {
            tags$ul(
              class = "small mb-0",
              lapply(avoid_ids, function(x) tags$li(x))
            )
          } else if (n_avoid > 5) {
            tags$span(class = "small", paste0(" (", avoid_ids[1], ", ", avoid_ids[2], ", ...)"))
          }
        )
      )
    })

    # -------------------------------------------------------------------------
    # Design probes
    # -------------------------------------------------------------------------
    observeEvent(input$design, {
      # Check we have a selection
      n_desired <- sum(unlist(values$selection_types) == "desired")
      if (n_desired == 0) {
        showNotification("Please select at least one desired target first.",
                         type = "error")
        return()
      }

      # Show progress
      showNotification("Designing probes...", id = "designing", duration = NULL)

      tryCatch({
        # Get desired and avoid IDs
        desired_ids <- names(values$selection_types)[unlist(values$selection_types) == "desired"]
        avoid_ids <- names(values$selection_types)[unlist(values$selection_types) == "avoid"]

        # Create target selection
        selection <- create_target_selection(values$trna_data, desired_ids, avoid_ids)

        # Determine region
        region <- switch(input$target_region,
          "all" = NULL,
          "5prime" = c(1, 25),
          "3prime" = c(-25, -1),  # Will be converted to actual positions
          "middle" = c(25, 50),
          "divergent" = "auto"
        )

        # Design probes
        if (input$target_region == "divergent" || length(avoid_ids) > 0) {
          # Use selective design
          result <- design_probes_selective(
            selection,
            values$trna_data,
            min_length = input$probe_length[1],
            max_length = input$probe_length[2],
            top_n = input$n_probes * 2  # Get extra to filter by Tm
          )
          probes <- result$probes
        } else {
          # Simple design for single target
          probes <- design_probes(
            desired_ids[1],
            values$trna_data,
            similarity_data = values$similarity_data,
            target_mode = "single",
            min_length = input$probe_length[1],
            max_length = input$probe_length[2],
            region = if (is.null(region)) NULL else region,
            top_n = input$n_probes * 2  # Get extra to filter by Tm
          )
        }

        # Filter by Tm range if we have results
        if (!is.null(probes) && nrow(probes) > 0 && "tm_nn" %in% names(probes)) {
          probes <- probes[probes$tm_nn >= input$tm_range[1] &
                          probes$tm_nn <= input$tm_range[2], ]
          # Limit to requested number
          if (nrow(probes) > input$n_probes) {
            probes <- probes[1:input$n_probes, ]
          }
        }

        # Store results
        local$probes <- probes

        removeNotification("designing")

        if (is.null(probes) || nrow(probes) == 0) {
          showNotification("No probes found matching the criteria. Try adjusting parameters.",
                           type = "warning", duration = 5)
        } else {
          showNotification(paste0("Designed ", nrow(probes), " probes!"),
                           type = "message", duration = 3)
        }

      }, error = function(e) {
        removeNotification("designing")
        showNotification(paste0("Error: ", e$message), type = "error", duration = 10)
      })
    })

    # -------------------------------------------------------------------------
    # Return probes
    # -------------------------------------------------------------------------
    return(list(
      probes = reactive(local$probes)
    ))

  })
}
