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
    design_params = NULL        # Design parameters from step 3
  )

  # ===========================================================================
  # Load data when organism changes
  # ===========================================================================

  observeEvent(input$organism, {
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

    compartment_label <- if (compartment == "nuclear") "nuclear/cytoplasmic" else "mitochondrial"
    showNotification(
      paste0("Showing ", nrow(values$trna_data), " ", compartment_label, " tRNAs"),
      type = "message",
      duration = 3
    )
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

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
    switch(values$wizard_step,
      wizardStep1UI("step1"),
      wizardStep2UI("step2"),
      wizardStep3UI("step3"),
      wizardStep4UI("step4")
    )
  })

  # ===========================================================================
  # Wizard navigation buttons
  # ===========================================================================

  output$wizard_nav <- renderUI({
    step <- values$wizard_step

    back_btn <- if (step > 1) {
      actionButton("wizard_back", "Back", class = "btn-secondary")
    } else {
      tags$div()  # Empty placeholder
    }

    next_btn <- if (step < 4) {
      # Disable next if conditions not met
      disabled <- switch(step,
        is.null(values$wizard_goal),  # Step 1: need goal
        length(values$wizard_selection$ids) == 0,  # Step 2: need selection
        FALSE  # Step 3: can always proceed
      )

      btn_label <- if (step == 3) "Design Probes" else "Next"
      btn <- actionButton("wizard_next", btn_label,
                          class = "btn-primary",
                          disabled = if (disabled) "disabled" else NULL)
      btn
    } else {
      # Step 4: Start Over button
      actionButton("wizard_reset", "Start Over", class = "btn-outline-primary")
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
    # This is the key metric - we want this to be low
    max_similarity <- divergence$max_cross_identity %||% 100

    # Determine status based on both conservation within targets AND specificity from non-targets
    # Lower max_similarity = easier to design specific probes
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
      return(list(status = "challenging", label = "Challenging",
                  message = "Some non-targets are similar; careful region selection needed"))
    } else {
      return(list(status = "not-feasible", label = "Not Feasible",
                  message = "Targets too similar to non-targets for specific detection"))
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

      values$wizard_probes <- result$probes
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
