# ui.R
# Compass Shiny app UI - Wizard-based probe design flow

ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#0072B2"  # Okabe-Ito blue
  ),

  # Include CSS
  tags$head(
    get_visualization_css(),
    tags$style(HTML("
      /* Wizard step indicator */
      .wizard-steps {
        display: flex;
        justify-content: center;
        margin: 20px 0 30px 0;
        padding: 0;
        list-style: none;
      }
      .wizard-step {
        display: flex;
        align-items: center;
        color: #999;
      }
      .wizard-step-number {
        width: 36px;
        height: 36px;
        border-radius: 50%;
        background-color: #e9ecef;
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: bold;
        margin-right: 8px;
      }
      .wizard-step-label {
        font-size: 0.9em;
      }
      .wizard-step.active .wizard-step-number {
        background-color: #0072B2;
        color: white;
      }
      .wizard-step.active .wizard-step-label {
        color: #0072B2;
        font-weight: bold;
      }
      .wizard-step.completed .wizard-step-number {
        background-color: #009E73;
        color: white;
      }
      .wizard-step.completed .wizard-step-label {
        color: #009E73;
      }
      .wizard-step-connector {
        width: 60px;
        height: 2px;
        background-color: #e9ecef;
        margin: 0 15px;
      }
      .wizard-step.completed + .wizard-step-connector,
      .wizard-step-connector.completed {
        background-color: #009E73;
      }

      /* Wizard content area */
      .wizard-content {
        min-height: 500px;
        padding: 20px;
      }

      /* Wizard navigation buttons */
      .wizard-nav {
        display: flex;
        justify-content: space-between;
        padding: 20px;
        border-top: 1px solid #dee2e6;
        margin-top: 20px;
      }

      /* Goal cards */
      .goal-card {
        cursor: pointer;
        transition: all 0.2s ease;
        border: 2px solid transparent;
        height: 100%;
      }
      .goal-card:hover {
        border-color: #0072B2;
        transform: translateY(-2px);
      }
      .goal-card.selected {
        border-color: #0072B2;
        background-color: #f0f7fb;
      }
      .goal-card .card-title {
        color: #0072B2;
      }

      /* Selection summary badge */
      .selection-badge {
        display: inline-block;
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.85em;
        margin-right: 10px;
      }
      .selection-badge.targets {
        background-color: #009E73;
        color: white;
      }

      /* Feasibility indicators */
      .feasibility-excellent { color: #009E73; }
      .feasibility-good { color: #56B4E9; }
      .feasibility-moderate { color: #E69F00; }
      .feasibility-challenging { color: #D55E00; }
      .feasibility-not-feasible { color: #CC79A7; }
    "))
  ),

  # Header
  tags$div(
    class = "container-fluid",
    style = "padding: 15px 20px; background-color: #0072B2; color: white;",
    tags$div(
      class = "d-flex justify-content-between align-items-center",
      tags$div(
        tags$h4("COMPASS", style = "margin: 0; font-weight: bold;"),
        tags$small("tRNA Northern Probe Designer")
      ),
      selectInput(
        "organism",
        label = NULL,
        choices = ORGANISMS,
        selected = "human",
        width = "150px"
      )
    )
  ),

  # Wizard steps indicator
  uiOutput("wizard_steps_ui"),

  # Main wizard content area
  tags$div(
    class = "container",
    tags$div(
      class = "wizard-content",
      uiOutput("wizard_content")
    ),

    # Navigation buttons
    uiOutput("wizard_nav")
  )
)
