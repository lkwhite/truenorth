# server.R
# TRUENORTH Shiny app server logic - Wizard-based flow

server <- function(input, output, session) {

  # ===========================================================================
  # Reactive values for app state
  # ===========================================================================

  values <- reactiveValues(
    # Data
    trna_data = NULL,           # Loaded tRNA data frame with anticodon positions
    similarity_data = NULL,     # Pre-computed similarity matrix

    # Wizard state
    wizard_step = 1,            # Current step (1-4)
    wizard_goal = NULL,         # "specific", "isoacceptor", "amino_acid", "consensus"

    # Selection state
    wizard_selection = list(
      type = NULL,              # Matches goal type
      ids = character(),        # Selected tRNA IDs
      amino_acid = NULL,        # For amino_acid goal
      anticodon = NULL          # For isoacceptor goal
    ),
    wizard_selection_obj = NULL,  # Selection object from create_target_selection

    # Analysis results
    wizard_feasibility = NULL,  # Feasibility analysis results
    wizard_probes = NULL,       # Designed probes
    design_params = NULL,       # Design parameters from step 3

    # Coverage tracking
    coverage_data = NULL,       # Coverage matrix and related data
    coverage_estimate = NULL,   # Estimated probes needed for coverage levels
    selected_probes = integer() # User-selected probe ranks for export
  )

  # ===========================================================================
  # Load data when organism changes
  # ===========================================================================

  observeEvent(input$organism, {
    req(input$organism)  # Don't run if organism is NULL
    showNotification("Loading tRNA data...", id = "loading", duration = NULL)

    # Load ALL tRNA data (both compartments) for off-target checking
    all_data <- load_trna_data(input$organism)
    values$trna_data_all <- all_data
    values$similarity_data <- load_sim_data(input$organism)

    # Filter by current compartment selection
    compartment <- input$compartment %||% "nuclear"
    values$trna_data <- all_data[all_data$compartment == compartment, ]

    # Reset wizard state when organism changes
    values$wizard_step <- 1
    values$wizard_goal <- NULL
    values$wizard_selection <- list(type = NULL, ids = character(), amino_acid = NULL, anticodon = NULL)
    values$wizard_selection_obj <- NULL
    values$wizard_feasibility <- NULL
    values$wizard_probes <- NULL
    values$design_params <- NULL
    values$coverage_data <- NULL
    values$coverage_estimate <- NULL
    values$selected_probes <- integer()

    removeNotification("loading")

    showNotification(
      paste0("Loaded ", nrow(values$trna_data), " ", compartment, " tRNAs for ",
             names(ORGANISMS)[ORGANISMS == input$organism]),
      type = "message",
      duration = 3
    )
  }, ignoreNULL = FALSE)

  # ===========================================================================
  # Filter by compartment when it changes
  # ===========================================================================

  observeEvent(input$compartment, {
    req(values$trna_data_all)

    compartment <- input$compartment
    values$trna_data <- values$trna_data_all[values$trna_data_all$compartment == compartment, ]

    # Reset wizard state when compartment changes
    values$wizard_step <- 1
    values$wizard_goal <- NULL
    values$wizard_selection <- list(type = NULL, ids = character(), amino_acid = NULL, anticodon = NULL)
    values$wizard_selection_obj <- NULL
    values$wizard_feasibility <- NULL
    values$wizard_probes <- NULL
    values$design_params <- NULL
    values$coverage_data <- NULL
    values$coverage_estimate <- NULL
    values$selected_probes <- integer()

    compartment_label <- if (compartment == "nuclear") "nuclear/cytoplasmic" else "mitochondrial"
    showNotification(
      paste0("Showing ", nrow(values$trna_data), " ", compartment_label, " tRNAs"),
      type = "message",
      duration = 3
    )
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # ===========================================================================
  # Organism display in header (shows after Step 1)
  # ===========================================================================

  output$organism_display <- renderUI({
    org <- input$organism %||% "human"
    comp <- input$compartment %||% "nuclear"

    org_label <- names(ORGANISMS)[ORGANISMS == org]
    comp_label <- if (comp == "nuclear") "Nuclear" else "Mito"

    if (values$wizard_step > 1) {
      # After step 1, show as static badges
      tags$div(
        class = "d-flex align-items-center gap-2",
        tags$span(
          class = "badge",
          style = "background-color: rgba(255,255,255,0.2); font-size: 0.9em;",
          org_label
        ),
        tags$span(
          class = "badge",
          style = "background-color: rgba(255,255,255,0.2); font-size: 0.9em;",
          comp_label
        )
      )
    }
    # On step 1, don't show anything in header (selectors are in the main content)
  })

  # ===========================================================================
  # Wizard steps indicator
  # ===========================================================================

  output$wizard_steps_ui <- renderUI({
    step <- values$wizard_step
    steps <- list(
      list(num = 1, label = "Goal"),
      list(num = 2, label = "Selection"),
      list(num = 3, label = "Feasibility"),
      list(num = 4, label = "Design")
    )

    step_elements <- lapply(seq_along(steps), function(i) {
      s <- steps[[i]]
      class <- "wizard-step"
      if (i < step) class <- paste(class, "completed")
      if (i == step) class <- paste(class, "active")

      connector <- NULL
      if (i < length(steps)) {
        conn_class <- if (i < step) "wizard-step-connector completed" else "wizard-step-connector"
        connector <- tags$div(class = conn_class)
      }

      tagList(
        tags$div(
          class = class,
          tags$div(class = "wizard-step-number",
                   if (i < step) icon("check") else s$num),
          tags$div(class = "wizard-step-label", s$label)
        ),
        connector
      )
    })

    tags$div(
      class = "container",
      tags$ul(class = "wizard-steps", step_elements)
    )
  })

  # ===========================================================================
  # Wizard content - renders current step
  # ===========================================================================

  output$wizard_content <- renderUI({
    if (values$wizard_step == 1) {
      # Step 1: Show organism/compartment selectors centered, then goal cards
      tagList(
        # Organism/compartment selection - centered
        tags$div(
          class = "text-center mb-4",
          tags$div(
            class = "d-inline-flex gap-3 align-items-center",
            selectInput(
              "organism",
              label = NULL,
              choices = ORGANISMS,
              selected = input$organism %||% "human",
              width = "150px"
            ),
            selectInput(
              "compartment",
              label = NULL,
              choices = c("Nuclear/Cytoplasmic" = "nuclear",
                          "Mitochondrial" = "mitochondrial"),
              selected = input$compartment %||% "nuclear",
              width = "180px"
            )
          )
        ),
        wizardStep1UI("step1")
      )
    } else {
      switch(values$wizard_step,
        NULL,  # Step 1 handled above
        wizardStep2UI("step2"),
        wizardStep3UI("step3"),
        wizardStep4UI("step4")
      )
    }
  })

  # ===========================================================================
  # Wizard navigation buttons
  # ===========================================================================

  output$wizard_nav <- renderUI({
    step <- values$wizard_step
    goal <- values$wizard_goal

    back_btn <- if (step > 1) {
      actionButton("wizard_back", "Back", class = "btn-secondary")
    } else {
      tags$div()  # Empty placeholder
    }

    # For Step 2 with isoacceptor/amino_acid goals, clicking auto-advances - no Next needed
    # For "specific" goal, user needs Next button to proceed after selecting multiple tRNAs
    hide_next_on_step2 <- step == 2 && goal %in% c("isoacceptor", "amino_acid", "consensus")

    next_btn <- if (step < 4 && !hide_next_on_step2) {
      # Disable next if conditions not met
      disabled <- switch(step,
        is.null(goal),  # Step 1: need goal
        length(values$wizard_selection$ids) == 0,  # Step 2: need selection
        FALSE  # Step 3: can always proceed
      )

      btn_label <- if (step == 3) "Design Probes" else "Next"
      btn <- actionButton("wizard_next", btn_label,
                          class = "btn-primary",
                          disabled = if (disabled) "disabled" else NULL)
      btn
    } else if (step == 4) {
      # Step 4: Start Over button
      actionButton("wizard_reset", "Start Over", class = "btn-outline-primary")
    } else {
      tags$div()  # Empty placeholder for Step 2 with auto-advance goals
    }

    tags$div(
      class = "wizard-nav",
      back_btn,
      next_btn
    )
  })

  # ===========================================================================
  # Navigation handlers
  # ===========================================================================

  observeEvent(input$wizard_back, {
    if (values$wizard_step > 1) {
      values$wizard_step <- values$wizard_step - 1
    }
  })

  observeEvent(input$wizard_next, {
    if (values$wizard_step < 4) {
      # Validate before advancing
      if (values$wizard_step == 1 && is.null(values$wizard_goal)) {
        showNotification("Please select a goal first", type = "warning")
        return()
      }
      if (values$wizard_step == 2 && length(values$wizard_selection$ids) == 0) {
        showNotification("Please make a selection first", type = "warning")
        return()
      }

      # Run feasibility analysis when moving to step 3
      if (values$wizard_step == 2) {
        run_feasibility_analysis()
      }

      # Run probe design when moving to step 4
      if (values$wizard_step == 3) {
        run_probe_design()
      }

      values$wizard_step <- values$wizard_step + 1
    }
  })

  observeEvent(input$wizard_reset, {
    values$wizard_step <- 1
    values$wizard_goal <- NULL
    values$wizard_selection <- list(type = NULL, ids = character(), amino_acid = NULL, anticodon = NULL)
    values$wizard_selection_obj <- NULL
    values$wizard_feasibility <- NULL
    values$wizard_probes <- NULL
    values$design_params <- NULL
    values$coverage_data <- NULL
    values$coverage_estimate <- NULL
    values$selected_probes <- integer()
  })

  # Reset wizard when logo is clicked
  observeEvent(input$logo_clicked, {
    values$wizard_step <- 1
    values$wizard_goal <- NULL
    values$wizard_selection <- list(type = NULL, ids = character(), amino_acid = NULL, anticodon = NULL)
    values$wizard_selection_obj <- NULL
    values$wizard_feasibility <- NULL
    values$wizard_probes <- NULL
    values$design_params <- NULL
    values$coverage_data <- NULL
    values$coverage_estimate <- NULL
    values$selected_probes <- integer()
  })

  # ===========================================================================
  # Auto-advance for isoacceptor/amino_acid selections
  # ===========================================================================

  # Track previous selection to detect changes
  previous_selection <- reactiveVal(character())

  observe({
    # Only trigger on Step 2 for isoacceptor, amino_acid, or consensus goals
    req(values$wizard_step == 2)
    req(values$wizard_goal %in% c("isoacceptor", "amino_acid", "consensus"))
    req(length(values$wizard_selection$ids) > 0)

    # Check if selection actually changed (avoid re-triggering on back navigation)
    current_ids <- values$wizard_selection$ids
    prev_ids <- previous_selection()

    if (!identical(sort(current_ids), sort(prev_ids))) {
      # Update previous selection tracker
      previous_selection(current_ids)

      # Run feasibility analysis and advance to Step 3
      run_feasibility_analysis()
      values$wizard_step <- 3
    }
  })

  # Reset previous selection when going back to Step 1 or changing organism
  observeEvent(values$wizard_step, {
    if (values$wizard_step == 1) {
      previous_selection(character())
    }
  })

  # ===========================================================================
  # Feasibility analysis
  # ===========================================================================

  run_feasibility_analysis <- function() {
    req(values$trna_data)
    req(length(values$wizard_selection$ids) > 0)

    showNotification("Analyzing feasibility...", id = "analyzing", duration = NULL)

    tryCatch({
      # Create selection object using the backend function
      selection <- create_target_selection(
        trna_df = values$trna_data,
        desired_ids = values$wizard_selection$ids
        # avoid_ids defaults to everything not in desired
      )

      # Get conservation analysis within desired group
      conservation <- analyze_group_conservation(
        selection$desired,
        values$similarity_data
      )

      # Get divergence from non-targets
      divergence <- analyze_group_divergence(
        selection,
        values$similarity_data
      )

      # Find best regions (can be slow, wrap in tryCatch)
      regions <- tryCatch({
        find_selective_regions(selection, min_length = 18, max_length = 22)
      }, error = function(e) {
        data.frame()  # Return empty if analysis fails
      })

      # Store selection for later use
      values$wizard_selection_obj <- selection

      values$wizard_feasibility <- list(
        n_targets = selection$n_desired,
        n_nontargets = selection$n_avoid,
        conservation = conservation,
        divergence = divergence,
        regions = regions,
        status = determine_feasibility_status(conservation, divergence)
      )

      removeNotification("analyzing")
    }, error = function(e) {
      removeNotification("analyzing")
      showNotification(paste("Analysis error:", e$message), type = "error")
    })
  }

  # Determine overall feasibility status
  determine_feasibility_status <- function(conservation, divergence) {
    # Handle conservation - use mean_identity from the analysis
    cons_pct <- conservation$mean_identity %||% 0

    # Handle divergence - check for no avoid targets
    if (!is.null(divergence$n_avoid) && divergence$n_avoid == 0) {
      return(list(status = "excellent", label = "Excellent",
                  message = "No cross-reactivity concerns - all tRNAs are targets"))
    }

    # max_cross_identity = how similar the CLOSEST non-target is (lower = more specific)
    # mean_cross_identity = average similarity to all non-targets
    max_similarity <- divergence$max_cross_identity %||% 100
    mean_similarity <- divergence$mean_cross_identity %||% 50

    # Use both max and mean to determine feasibility
    # High max but low mean = one close relative (common for same amino acid family) - still workable
    # High max AND high mean = genuinely hard to distinguish

    if (max_similarity <= 70) {
      # Non-targets are quite different - good specificity potential
      if (cons_pct >= 90) {
        return(list(status = "excellent", label = "Excellent",
                    message = "High target conservation and good specificity"))
      } else if (cons_pct >= 70) {
        return(list(status = "good", label = "Good",
                    message = "Good specificity; some target variation may require multiple probes"))
      } else {
        return(list(status = "moderate", label = "Moderate",
                    message = "Good specificity but targets vary; may need multiple probes"))
      }
    } else if (max_similarity <= 80) {
      # Some similar non-targets but still workable
      if (cons_pct >= 85) {
        return(list(status = "good", label = "Good",
                    message = "Probe design should work well"))
      } else {
        return(list(status = "moderate", label = "Moderate",
                    message = "May need careful region selection"))
      }
    } else if (max_similarity <= 90) {
      return(list(status = "moderate", label = "Moderate",
                  message = "Some related tRNAs are similar; probe design will find distinguishing regions"))
    } else if (max_similarity <= 97) {
      # High max similarity - but check mean to see if it's just one close relative
      if (mean_similarity <= 50) {
        # Low mean = just a few close relatives (e.g., same amino acid family)
        return(list(status = "moderate", label = "Moderate",
                    message = "Close relatives exist (likely same amino acid); probe design can usually distinguish"))
      } else {
        return(list(status = "challenging", label = "Challenging",
                    message = "Several similar non-targets; careful region selection needed"))
      }
    } else {
      # >97% similar - nearly identical sequences
      return(list(status = "not-feasible", label = "Very Challenging",
                  message = "Very high similarity to some non-targets; distinguishing may be difficult"))
    }
  }

  # ===========================================================================
  # Probe design
  # ===========================================================================

  run_probe_design <- function() {
    req(values$trna_data)
    req(length(values$wizard_selection$ids) > 0)

    showNotification("Designing probes...", id = "designing", duration = NULL)

    tryCatch({
      # Get or create selection object
      selection <- values$wizard_selection_obj
      if (is.null(selection)) {
        selection <- create_target_selection(
          trna_df = values$trna_data,
          desired_ids = values$wizard_selection$ids
        )
      }

      # Get design parameters from step 3 (with defaults)
      params <- values$design_params %||% list(
        probe_length = 20,
        region_pref = "any",
        avoid_anticodon = TRUE
      )

      # Design probes using the backend function
      result <- design_probes_selective(
        selection = selection,
        trna_df = values$trna_data,
        min_length = params$probe_length - 2,
        max_length = params$probe_length + 2,
        min_conservation = 80,
        min_divergence = 20,
        top_n = 20,
        region_pref = params$region_pref,
        avoid_anticodon = params$avoid_anticodon
      )

      # Build coverage matrix for the designed probes
      targets_df <- values$trna_data[values$trna_data$id %in% values$wizard_selection$ids, ]

      if (nrow(result$probes) > 0 && nrow(targets_df) > 0) {
        coverage_result <- build_coverage_matrix(
          probes = result$probes,
          targets_df = targets_df,
          max_mismatches = 3
        )

        probes_final <- coverage_result$probes
        coverage_matrix <- coverage_result$matrix

        # For isoacceptor and amino_acid goals, re-rank probes for optimal coverage
        # This ensures rank 1, 2, 3... are complementary probes, not variations
        if (values$wizard_goal %in% c("isoacceptor", "amino_acid")) {
          probes_final <- rerank_probes_for_coverage(
            probes = probes_final,
            coverage_matrix = coverage_matrix,
            target_ids = targets_df$id
          )

          # Rebuild coverage matrix with new probe order
          coverage_result <- build_coverage_matrix(
            probes = probes_final,
            targets_df = targets_df,
            max_mismatches = 3
          )
          coverage_matrix <- coverage_result$matrix
          probes_final <- coverage_result$probes
        }

        # Store final probes
        values$wizard_probes <- probes_final

        # Store coverage data for UI
        values$coverage_data <- list(
          matrix = coverage_matrix,
          mismatch_matrix = coverage_result$mismatch_matrix,
          targets_covered_by = coverage_result$targets_covered_by,
          n_targets = coverage_result$n_targets,
          target_ids = targets_df$id
        )

        # Estimate probes needed for various coverage levels
        values$coverage_estimate <- estimate_probes_needed(
          coverage_matrix = coverage_matrix,
          target_ids = targets_df$id
        )
      } else {
        values$wizard_probes <- result$probes
        values$coverage_data <- NULL
        values$coverage_estimate <- NULL
      }

      removeNotification("designing")

    }, error = function(e) {
      removeNotification("designing")
      showNotification(paste("Design error:", e$message), type = "error")
    })
  }

  # ===========================================================================
  # Module servers
  # ===========================================================================

  wizardStep1Server("step1", values)
  wizardStep2Server("step2", values)
  wizardStep3Server("step3", values)
  wizardStep4Server("step4", values)

}
