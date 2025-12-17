# ui.R
# Compass Shiny app UI with sidebar navigation

ui <- page_navbar(
  title = "COMPASS",
  id = "navbar",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#0072B2"  # Okabe-Ito blue

),

  # Header with organism selector
  header = tags$div(
    class = "container-fluid",
    style = "padding: 10px 20px; background-color: #f8f9fa; border-bottom: 1px solid #dee2e6;",
    tags$div(
      class = "d-flex justify-content-between align-items-center",
      tags$span(
        style = "font-size: 0.9em; color: #666;",
        "tRNA Northern Probe Designer"
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

  # Include visualization CSS
  tags$head(
    get_visualization_css(),
    tags$style(HTML("
      .nav-link.active {
        font-weight: bold;
      }
      .selection-summary {
        padding: 10px;
        background-color: #f8f9fa;
        border-radius: 4px;
        margin-top: 15px;
      }
      .btn-desired {
        background-color: #009E73;
        border-color: #009E73;
        color: white;
      }
      .btn-desired:hover {
        background-color: #007d5c;
        border-color: #007d5c;
        color: white;
      }
      .btn-avoid {
        background-color: #D55E00;
        border-color: #D55E00;
        color: white;
      }
      .btn-avoid:hover {
        background-color: #b34d00;
        border-color: #b34d00;
        color: white;
      }
    "))
  ),

  # ==========================================================================
  # Browse tRNAs tab
  # ==========================================================================
  nav_panel(
    title = "Browse tRNAs",
    icon = icon("dna"),
    browseUI("browse")
  ),

  # ==========================================================================
  # Design Probes tab
  # ==========================================================================
  nav_panel(
    title = "Design Probes",
    icon = icon("flask"),
    designUI("design")
  ),

  # ==========================================================================
  # Results tab
  # ==========================================================================
  nav_panel(
    title = "Results",
    icon = icon("table"),
    resultsUI("results")
  ),

  # ==========================================================================
  # Sidebar with selection summary
  # ==========================================================================
  nav_spacer(),
  nav_item(
    tags$div(
      style = "padding: 10px; min-width: 150px;",
      tags$strong("Selection"),
      tags$div(
        id = "selection-summary",
        class = "selection-summary",
        tags$div(
          id = "desired-count",
          style = paste0("color: ", COLORS$desired, ";"),
          "0 desired"
        ),
        tags$div(
          id = "avoid-count",
          style = paste0("color: ", COLORS$avoid, ";"),
          "0 avoid"
        )
      )
    )
  )
)
