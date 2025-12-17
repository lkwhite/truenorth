# wizard_step1.R
# Step 1: Goal Selection - "What do you want to detect?"

wizardStep1UI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$div(
      class = "text-center mb-4",
      tags$h3("What do you want to detect?"),
      tags$p(class = "text-muted",
             "Choose your detection goal to help guide probe design")
    ),

    fluidRow(
      # Specific isodecoder
      column(
        width = 4,
        tags$div(
          id = ns("card_specific"),
          class = "card goal-card",
          onclick = sprintf("Shiny.setInputValue('%s', 'specific', {priority: 'event'})",
                            ns("goal_click")),
          tags$div(
            class = "card-body text-center",
            tags$div(
              style = "font-size: 2.5em; margin-bottom: 10px;",
              icon("crosshairs")
            ),
            tags$h5(class = "card-title", "Specific Isodecoder(s)"),
            tags$p(class = "card-text text-muted",
                   "Target specific tRNA gene variants while avoiding others"),
            tags$small(class = "text-muted",
                       "Best for: Variant-specific detection, distinguishing similar genes")
          )
        )
      ),

      # Isoacceptor family
      column(
        width = 4,
        tags$div(
          id = ns("card_isoacceptor"),
          class = "card goal-card",
          onclick = sprintf("Shiny.setInputValue('%s', 'isoacceptor', {priority: 'event'})",
                            ns("goal_click")),
          tags$div(
            class = "card-body text-center",
            tags$div(
              style = "font-size: 2.5em; margin-bottom: 10px;",
              icon("layer-group")
            ),
            tags$h5(class = "card-title", "Isoacceptor Family"),
            tags$p(class = "card-text text-muted",
                   "Target all tRNAs sharing an anticodon"),
            tags$small(class = "text-muted",
                       "Best for: Codon-specific detection (may need 1+ probes)")
          )
        )
      ),

      # All for amino acid
      column(
        width = 4,
        tags$div(
          id = ns("card_amino_acid"),
          class = "card goal-card",
          onclick = sprintf("Shiny.setInputValue('%s', 'amino_acid', {priority: 'event'})",
                            ns("goal_click")),
          tags$div(
            class = "card-body text-center",
            tags$div(
              style = "font-size: 2.5em; margin-bottom: 10px;",
              icon("shapes")
            ),
            tags$h5(class = "card-title", "Amino Acid Pool"),
            tags$p(class = "card-text text-muted",
                   "Target all tRNAs for an amino acid across anticodons"),
            tags$small(class = "text-muted",
                       "Best for: Total pool measurement (often needs multiple probes)")
          )
        )
      )
    ),

    # Selection feedback
    uiOutput(ns("selection_feedback"))
  )
}

wizardStep1Server <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Handle goal selection
    observeEvent(input$goal_click, {
      values$wizard_goal <- input$goal_click

      # Reset selection when goal changes
      values$wizard_selection <- list(
        type = input$goal_click,
        ids = character(),
        amino_acid = NULL,
        anticodon = NULL
      )

      # Update card styling via JavaScript
      session$sendCustomMessage("updateGoalCards", list(
        selected = input$goal_click,
        ns = ns("")
      ))
    })

    # Show selection feedback
    output$selection_feedback <- renderUI({
      req(values$wizard_goal)

      goal_labels <- list(
        specific = "Specific Isodecoder(s)",
        isoacceptor = "Isoacceptor Family",
        amino_acid = "Amino Acid Pool"
      )

      goal_descriptions <- list(
        specific = "You'll select which tRNA genes to target. We'll analyze feasibility and design probes to hit your targets while avoiding others.",
        isoacceptor = "You'll select an anticodon family. We'll analyze if one probe can hit all members, or if you need multiple probes.",
        amino_acid = "You'll select an amino acid. We'll determine how many probes are needed to cover all its tRNA isoacceptors."
      )

      tags$div(
        class = "alert alert-info mt-4",
        tags$strong(goal_labels[[values$wizard_goal]], " selected. "),
        goal_descriptions[[values$wizard_goal]],
        tags$br(),
        tags$small("Click Next to continue to target selection.")
      )
    })

  })
}
