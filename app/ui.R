# ui.R
# TRUENORTH Shiny app UI - Wizard-based probe design flow

ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#0072B2"  # Okabe-Ito blue
  ),

  # Include CSS and JS
  tags$head(
    get_visualization_css(),
    tooltip_js(),  # Initialize Bootstrap tooltips
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
    ")),
    # JavaScript for goal card selection
    tags$script(HTML("
      $(document).ready(function() {
        Shiny.addCustomMessageHandler('updateGoalCards', function(data) {
          // Remove selected class from all goal cards
          $('.goal-card').removeClass('selected');
          // Add selected class to chosen card
          var selectedCard = document.getElementById(data.ns + 'card_' + data.selected);
          if (selectedCard) {
            selectedCard.classList.add('selected');
          }
        });
      });
    "))
  ),

  # Header
  tags$div(
    class = "container-fluid",
    style = "padding: 15px 20px; background-color: #0072B2; color: white;",
    tags$div(
      class = "d-flex justify-content-between align-items-center",
      # Clickable logo + title to reset wizard
      tags$a(
        id = "logo_reset",
        href = "#",
        class = "text-decoration-none",
        style = "display: flex; align-items: center; gap: 12px; color: white; cursor: pointer;",
        onclick = "Shiny.setInputValue('logo_clicked', Math.random()); return false;",
        # Compass SVG logo - black/white with red needle pointing north
        HTML('
          <svg width="40" height="40" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">
            <!-- Outer circle -->
            <circle cx="50" cy="50" r="45" fill="none" stroke="white" stroke-width="3"/>
            <!-- Inner circle -->
            <circle cx="50" cy="50" r="35" fill="none" stroke="white" stroke-width="1.5"/>
            <!-- Cardinal directions -->
            <text x="50" y="18" text-anchor="middle" fill="white" font-size="10" font-weight="bold">N</text>
            <text x="50" y="92" text-anchor="middle" fill="white" font-size="10" font-weight="bold">S</text>
            <text x="8" y="54" text-anchor="middle" fill="white" font-size="10" font-weight="bold">W</text>
            <text x="92" y="54" text-anchor="middle" fill="white" font-size="10" font-weight="bold">E</text>
            <!-- Compass needle - red north, white south -->
            <polygon points="50,22 45,50 50,45 55,50" fill="#D55E00"/>
            <polygon points="50,78 45,50 50,55 55,50" fill="white"/>
            <!-- Center dot -->
            <circle cx="50" cy="50" r="4" fill="white"/>
          </svg>
        '),
        tags$div(
          tags$h4("TRUENORTH", style = "margin: 0; font-weight: bold;"),
          tags$small("tRNA-Resolved Utility for Evaluating Northern Oligo Reference Targets via Hybridization",
                     style = "opacity: 0.9;")
        )
      ),
      # Show organism/compartment selection as static text (after Step 1)
      uiOutput("organism_display")
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
