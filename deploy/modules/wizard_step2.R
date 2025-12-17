# wizard_step2.R
# Step 2: Target Selection - varies based on goal

wizardStep2UI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(
      class = "text-center mb-4",
      tags$h3("Select your targets"),
      uiOutput(ns("goal_subtitle"))
    ),

    # Dynamic content based on goal
    uiOutput(ns("selection_ui")),

    # Selection summary
    uiOutput(ns("selection_summary"))
  )
}

wizardStep2Server <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Subtitle based on goal
    output$goal_subtitle <- renderUI({
      goal <- values$wizard_goal
      subtitle <- switch(goal,
        "specific" = "Click on tRNAs to select specific isodecoders",
        "isoacceptor" = "Select an anticodon family",
        "amino_acid" = "Select an amino acid to target all its tRNAs",
        "consensus" = "Review consensus probe targets",
        "Choose your targets"
      )
      tags$p(class = "text-muted", subtitle)
    })

    # Dynamic selection UI based on goal
    output$selection_ui <- renderUI({
      req(values$wizard_goal)
      req(values$trna_data)

      switch(values$wizard_goal,
        "specific" = render_specific_ui(ns, values$trna_data),
        "isoacceptor" = render_isoacceptor_ui(ns, values$trna_data),
        "amino_acid" = render_amino_acid_ui(ns, values$trna_data),
        "consensus" = render_consensus_ui(ns, values$trna_data)
      )
    })

    # Update anticodon filter when amino acid filter changes (for specific selection)
    observeEvent(input$filter_amino_acid, {
      req(values$trna_data)

      if (is.null(input$filter_amino_acid) || input$filter_amino_acid == "") {
        # Show all anticodons
        anticodons <- sort(unique(values$trna_data$anticodon))
      } else {
        # Filter anticodons by selected amino acid
        filtered <- values$trna_data[values$trna_data$amino_acid == input$filter_amino_acid, ]
        anticodons <- sort(unique(filtered$anticodon))
      }

      updateSelectInput(
        session,
        "filter_anticodon",
        choices = c("All" = "", anticodons),
        selected = ""
      )
    }, ignoreNULL = FALSE)

    # Update tRNA checkbox list when filters change
    observe({
      req(values$trna_data)
      req(values$wizard_goal == "specific")

      filtered_data <- values$trna_data

      # Apply amino acid filter
      if (!is.null(input$filter_amino_acid) && input$filter_amino_acid != "") {
        filtered_data <- filtered_data[filtered_data$amino_acid == input$filter_amino_acid, ]
      }

      # Apply anticodon filter
      if (!is.null(input$filter_anticodon) && input$filter_anticodon != "") {
        filtered_data <- filtered_data[filtered_data$anticodon == input$filter_anticodon, ]
      }

      # Update checkbox choices while preserving selections
      current_selection <- input$selected_trnas %||% character()

      updateCheckboxGroupInput(
        session,
        "selected_trnas",
        choices = setNames(filtered_data$id, format_trna_id(filtered_data$id)),
        selected = intersect(current_selection, filtered_data$id)
      )
    })

    # Handle specific isodecoder selection
    observeEvent(input$selected_trnas, {
      if (!is.null(values$wizard_goal) && values$wizard_goal == "specific") {
        values$wizard_selection$ids <- input$selected_trnas %||% character()
      }
    }, ignoreNULL = FALSE)

    # Handle isoacceptor family selection
    observeEvent(input$selected_anticodon, {
      if (!is.null(values$wizard_goal) && values$wizard_goal == "isoacceptor" &&
          !is.null(input$selected_anticodon) && !is.null(values$trna_data)) {
        values$wizard_selection$anticodon <- input$selected_anticodon
        # Get all tRNA IDs with this anticodon
        matching <- values$trna_data[values$trna_data$anticodon == input$selected_anticodon, ]
        values$wizard_selection$ids <- matching$id
      }
    }, ignoreNULL = FALSE)

    # Handle amino acid selection
    observeEvent(input$selected_amino_acid, {
      if (!is.null(values$wizard_goal) && values$wizard_goal == "amino_acid" &&
          !is.null(input$selected_amino_acid) && !is.null(values$trna_data)) {
        values$wizard_selection$amino_acid <- input$selected_amino_acid
        # Get all tRNA IDs for this amino acid
        matching <- values$trna_data[values$trna_data$amino_acid == input$selected_amino_acid, ]
        values$wizard_selection$ids <- matching$id
      }
    }, ignoreNULL = FALSE)

    # Handle consensus auto-selection
    observe({
      if (!is.null(values$wizard_goal) && values$wizard_goal == "consensus" &&
          !is.null(values$trna_data)) {
        # Auto-select all tRNAs for consensus analysis
        values$wizard_selection$ids <- values$trna_data$id
      }
    })

    # Selection summary
    output$selection_summary <- renderUI({
      n_selected <- length(values$wizard_selection$ids)

      if (n_selected == 0) {
        return(tags$div(
          class = "alert alert-warning mt-3",
          "No targets selected. Please make a selection to continue."
        ))
      }

      # Get summary stats
      selected_data <- values$trna_data[values$trna_data$id %in% values$wizard_selection$ids, ]
      n_anticodons <- length(unique(selected_data$anticodon))
      n_amino_acids <- length(unique(selected_data$amino_acid))

      tags$div(
        class = "alert alert-success mt-3",
        tags$strong(n_selected, " tRNAs selected"),
        tags$br(),
        tags$small(
          n_anticodons, " anticodon families, ",
          n_amino_acids, " amino acids"
        )
      )
    })
  })
}

# =============================================================================
# Helper functions for rendering different selection UIs
# =============================================================================

render_specific_ui <- function(ns, trna_data) {
  # Group by amino acid and anticodon for organized display
  aa_choices <- sort(unique(trna_data$amino_acid))

  tagList(
    fluidRow(
      column(
        width = 4,
        selectInput(
          ns("filter_amino_acid"),
          "Filter by amino acid:",
          choices = c("All" = "", aa_choices),
          selected = ""
        )
      ),
      column(
        width = 4,
        selectInput(
          ns("filter_anticodon"),
          "Filter by anticodon:",
          choices = c("All" = ""),
          selected = ""
        )
      )
    ),

    tags$div(
      style = "max-height: 400px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 4px; padding: 10px;",
      checkboxGroupInput(
        ns("selected_trnas"),
        label = NULL,
        choices = setNames(trna_data$id, format_trna_id(trna_data$id)),
        selected = character()
      )
    )
  )
}

render_isoacceptor_ui <- function(ns, trna_data) {
  # Get anticodon families with counts
  anticodon_summary <- trna_data %>%
    group_by(amino_acid, anticodon) %>%
    summarise(
      n_genes = n(),
      .groups = "drop"
    ) %>%
    arrange(amino_acid, anticodon)

  # Create card grid for anticodons
  cards <- lapply(seq_len(nrow(anticodon_summary)), function(i) {
    row <- anticodon_summary[i, ]
    card_id <- paste0("ac_card_", row$anticodon)

    tags$div(
      class = "col-lg-2 col-md-3 col-sm-4 col-6 mb-2",
      tags$div(
        id = ns(card_id),
        class = "card goal-card anticodon-card",
        style = "cursor: pointer;",
        onclick = sprintf(
          "document.querySelectorAll('.anticodon-card').forEach(c => c.classList.remove('selected')); this.classList.add('selected'); Shiny.setInputValue('%s', '%s', {priority: 'event'})",
          ns("selected_anticodon"), row$anticodon
        ),
        tags$div(
          class = "card-body text-center py-2",
          tags$h6(class = "card-title mb-0", row$anticodon),
          tags$small(class = "text-muted", row$amino_acid),
          tags$div(
            tags$span(class = "badge bg-secondary", style = "font-size: 0.7em;", paste(row$n_genes, "genes"))
          )
        )
      )
    )
  })

  tags$div(
    class = "row",
    cards
  )
}

render_amino_acid_ui <- function(ns, trna_data) {
  # Get amino acid summary
  aa_summary <- trna_data %>%
    group_by(amino_acid) %>%
    summarise(
      n_genes = n(),
      n_anticodons = n_distinct(anticodon),
      .groups = "drop"
    ) %>%
    arrange(amino_acid)

  # Single-letter amino acid codes
  aa_codes <- c(
    "Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Cys" = "C",
    "Gln" = "Q", "Glu" = "E", "Gly" = "G", "His" = "H", "Ile" = "I",
    "Leu" = "L", "Lys" = "K", "Met" = "M", "Phe" = "F", "Pro" = "P",
    "Ser" = "S", "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V",
    "SeC" = "U", "Sup" = "*", "iMet" = "iM", "Und" = "?"
  )

  cards <- lapply(seq_len(nrow(aa_summary)), function(i) {
    row <- aa_summary[i, ]
    code <- aa_codes[row$amino_acid]
    if (is.na(code)) code <- substr(row$amino_acid, 1, 1)

    tags$div(
      class = "col-lg-1 col-md-2 col-sm-3 col-4 mb-2",
      tags$div(
        id = ns(paste0("aa_card_", row$amino_acid)),
        class = "card goal-card amino-acid-card",
        style = "cursor: pointer;",
        onclick = sprintf(
          "document.querySelectorAll('.amino-acid-card').forEach(c => c.classList.remove('selected')); this.classList.add('selected'); Shiny.setInputValue('%s', '%s', {priority: 'event'})",
          ns("selected_amino_acid"), row$amino_acid
        ),
        tags$div(
          class = "card-body text-center py-1 px-1",
          tags$div(
            style = "font-size: 1.5em; font-weight: bold; color: #0072B2;",
            code
          ),
          tags$small(class = "d-block", style = "font-size: 0.75em;", row$amino_acid),
          tags$small(class = "text-muted", style = "font-size: 0.7em;",
                     paste(row$n_genes, "genes"))
        )
      )
    )
  })

  tags$div(
    class = "row justify-content-center",
    cards
  )
}

render_consensus_ui <- function(ns, trna_data) {
  n_total <- nrow(trna_data)
  n_anticodons <- length(unique(trna_data$anticodon))
  n_amino_acids <- length(unique(trna_data$amino_acid))

  tags$div(
    class = "text-center",
    tags$div(
      class = "card",
      style = "max-width: 500px; margin: 0 auto;",
      tags$div(
        class = "card-body",
        tags$div(
          style = "font-size: 3em; color: #0072B2;",
          icon("bullseye")
        ),
        tags$h4("Consensus Probe Analysis"),
        tags$p(class = "text-muted",
               "We'll analyze all tRNAs to find conserved regions that can serve as pan-tRNA probes."),
        tags$hr(),
        tags$div(
          class = "row text-center",
          tags$div(
            class = "col",
            tags$h3(n_total),
            tags$small(class = "text-muted", "Total tRNAs")
          ),
          tags$div(
            class = "col",
            tags$h3(n_anticodons),
            tags$small(class = "text-muted", "Anticodons")
          ),
          tags$div(
            class = "col",
            tags$h3(n_amino_acids),
            tags$small(class = "text-muted", "Amino Acids")
          )
        ),
        tags$hr(),
        tags$p(
          class = "mb-0",
          tags$small(
            "Click ", tags$strong("Next"), " to analyze consensus regions across all tRNAs."
          )
        )
      )
    )
  )
}
