# help_content.R
# In-app help content and tooltips for TRUENORTH
#
# This module provides contextual help content for each wizard step,
# terminology definitions, and reusable help UI components.

library(shiny)
library(htmltools)

# =============================================================================
# Terminology Definitions
# =============================================================================

TERMINOLOGY <- list(
  isoacceptor = list(
    term = "Isoacceptor",
    short = "tRNAs sharing the same anticodon",
    full = "All tRNAs that share the same anticodon sequence. They decode the same codon(s) and carry the same amino acid. For example, all tRNA-Ala-AGC molecules are isoacceptors."
  ),
  isodecoder = list(
    term = "Isodecoder",
    short = "Same anticodon, different sequence",
    full = "tRNAs with the same anticodon but different body sequences. Even though they decode the same codons, their sequence differences can be exploited for specific probe design. Isodecoders may not cross-hybridize to the same probe."
  ),
  gene_family = list(
    term = "Gene family",
    short = "Copies from same ancestral gene",
    full = "tRNA genes derived from the same ancestral gene through duplication. Family members typically have identical or nearly identical sequences. The naming convention is AA-anticodon-family-copy (e.g., Ala-AGC-1-1, Ala-AGC-1-2 are copies within family 1)."
  ),
  anticodon = list(
    term = "Anticodon",
    short = "3-nucleotide codon-binding sequence",
    full = "The three-nucleotide sequence in the tRNA that base-pairs with the mRNA codon during translation. Located in the anticodon loop around positions 34-36."
  ),
  tm = list(
    term = "Melting Temperature (Tm)",
    short = "Temperature at which 50% of probe is hybridized",
    full = "The temperature at which half of the probe molecules are bound to their target. Calculated using nearest-neighbor thermodynamics. Higher Tm indicates stronger binding."
  ),
  gc_content = list(
    term = "GC Content",
    short = "Percentage of G and C bases",
    full = "The percentage of guanine and cytosine bases in the probe sequence. Higher GC content generally increases melting temperature and binding stability. 40-60% is optimal for most applications."
  ),
  modification = list(
    term = "Modified Regions",
    short = "Heavily chemically modified positions",
    full = "tRNAs contain many post-transcriptional modifications. Key modified regions include the D-loop (dihydrouridine), anticodon loop (multiple modifications), and TΨC loop (pseudouridine). Probes targeting these regions may have unpredictable hybridization."
  )
)

# =============================================================================
# Tooltip Helper
# =============================================================================

#' Create a tooltip icon with term definition
#'
#' @param term_key Key from TERMINOLOGY list (e.g., "isoacceptor")
#' @param style Optional additional CSS style
#' @return Span element with tooltip
create_tooltip <- function(term_key, style = NULL) {
  if (!term_key %in% names(TERMINOLOGY)) {
    return(tags$span())
  }

  term <- TERMINOLOGY[[term_key]]

  tags$span(
    class = "help-tooltip",
    style = paste0("cursor: help; border-bottom: 1px dotted #666; ", style),
    `data-bs-toggle` = "tooltip",
    `data-bs-placement` = "top",
    `data-bs-html` = "true",
    title = paste0("<strong>", term$term, "</strong><br>", term$short),
    term$term
  )
}

#' Create an inline help icon (question mark)
#'
#' @param help_text Text to show on hover
#' @return Span element with help icon
help_icon <- function(help_text) {
  tags$span(
    class = "help-icon",
    style = "cursor: help; color: #0072B2; margin-left: 5px;",
    `data-bs-toggle` = "tooltip",
    `data-bs-placement` = "top",
    `data-bs-html` = "true",
    title = help_text,
    HTML("&#9432;")  # Information icon
  )
}

# =============================================================================
# Step-specific Help Content
# =============================================================================

#' Help content for Step 1: Goal Selection
help_step1_ui <- function() {
  tags$div(
    class = "help-section",
    style = "background-color: #f8f9fa; border-radius: 8px; padding: 15px; margin-bottom: 20px;",

    tags$div(
      style = "display: flex; align-items: center; margin-bottom: 10px;",
      tags$span(style = "font-size: 1.2em; margin-right: 8px;", HTML("&#128161;")),
      tags$strong("Choosing Your Detection Goal")
    ),

    tags$p(
      style = "margin-bottom: 10px; color: #555;",
      "Your goal determines how specific your probe needs to be:"
    ),

    tags$ul(
      style = "margin-bottom: 0; color: #555;",
      tags$li(
        tags$strong("Specific Isodecoder(s):"),
        " When you need to distinguish between tRNAs with the same anticodon but different sequences. Most challenging but highest specificity."
      ),
      tags$li(
        tags$strong("Isoacceptor Family:"),
        " When you want to detect all tRNAs with a particular anticodon (e.g., all Ala-AGC). Good balance of specificity and coverage."
      ),
      tags$li(
        tags$strong("Amino Acid Group:"),
        " When you want to detect all tRNAs for an amino acid regardless of anticodon. Broadest detection, easier probe design."
      )
    )
  )
}

#' Help content for Step 2: Target Selection
help_step2_ui <- function() {
  tags$div(
    class = "help-section",
    style = "background-color: #f8f9fa; border-radius: 8px; padding: 15px; margin-bottom: 20px;",

    tags$div(
      style = "display: flex; align-items: center; margin-bottom: 10px;",
      tags$span(style = "font-size: 1.2em; margin-right: 8px;", HTML("&#128161;")),
      tags$strong("Selecting Your Targets")
    ),

    tags$p(
      style = "margin-bottom: 10px; color: #555;",
      "Click on tiles to select the tRNAs you want your probe to detect. Multiple selections are allowed."
    ),

    tags$div(
      style = "font-size: 0.9em; color: #666;",
      tags$p(
        style = "margin-bottom: 5px;",
        tags$strong("Tip: "),
        "Identical sequences are grouped together. Selecting a group selects all members."
      ),
      tags$p(
        style = "margin-bottom: 0;",
        tags$strong("Note: "),
        "The highlighted ",
        tags$span(style = "background-color: #000; color: white; padding: 0 4px;", "NNN"),
        " region shows the anticodon position."
      )
    )
  )
}

#' Help content for Step 3: Feasibility Review
help_step3_ui <- function() {
  tags$div(
    class = "help-section",
    style = "background-color: #f8f9fa; border-radius: 8px; padding: 15px; margin-bottom: 20px;",

    tags$div(
      style = "display: flex; align-items: center; margin-bottom: 10px;",
      tags$span(style = "font-size: 1.2em; margin-right: 8px;", HTML("&#128161;")),
      tags$strong("Understanding Feasibility")
    ),

    tags$p(
      style = "margin-bottom: 10px; color: #555;",
      "This step analyzes whether your target selection allows for specific probe design."
    ),

    tags$ul(
      style = "margin-bottom: 10px; color: #555; font-size: 0.9em;",
      tags$li(
        tags$strong("Conservation:"),
        " How similar are your targets to each other? Higher conservation means one probe can detect them all."
      ),
      tags$li(
        tags$strong("Divergence:"),
        " How different are your targets from non-targets? Higher divergence means better specificity."
      ),
      tags$li(
        tags$strong("Closest Off-targets:"),
        " The non-target tRNAs most likely to cross-hybridize. Watch for high identity scores."
      )
    ),

    tags$div(
      style = "padding: 10px; background-color: #fff3cd; border-radius: 4px; font-size: 0.85em;",
      tags$strong("Challenging detection? "),
      "Consider narrowing your targets or accepting that some cross-hybridization may occur."
    )
  )
}

#' Help content for Step 4: Probe Design
help_step4_ui <- function() {
  tags$div(
    class = "help-section",
    style = "background-color: #f8f9fa; border-radius: 8px; padding: 15px; margin-bottom: 20px;",

    tags$div(
      style = "display: flex; align-items: center; margin-bottom: 10px;",
      tags$span(style = "font-size: 1.2em; margin-right: 8px;", HTML("&#128161;")),
      tags$strong("Selecting Your Probes")
    ),

    tags$p(
      style = "margin-bottom: 10px; color: #555;",
      "Probes are ranked by combined quality score. Select probes using the checkboxes to build your probe set."
    ),

    tags$ul(
      style = "margin-bottom: 10px; color: #555; font-size: 0.9em;",
      tags$li(
        tags$strong("Coverage tracker:"),
        " Shows how many targets your selected probes detect. Aim for 100% coverage."
      ),
      tags$li(
        tags$strong("Modification warnings:"),
        " ",
        tags$span(style = "color: #E69F00;", HTML("&#9888;")),
        " indicates overlap with heavily modified regions. Use with caution."
      ),
      tags$li(
        tags$strong("Binding diagrams:"),
        " Click a row to expand and see exactly how the probe binds to targets."
      )
    ),

    tags$div(
      style = "padding: 10px; background-color: #e7f3ff; border-radius: 4px; font-size: 0.85em;",
      tags$strong("Tip: "),
      "If one probe doesn't cover all targets, select additional probes to increase coverage. The coverage tracker updates in real-time."
    )
  )
}

# =============================================================================
# Collapsible Help Panel
# =============================================================================

#' Create a collapsible help panel
#'
#' @param id Unique ID for the panel
#' @param title Panel title
#' @param content HTML content to display
#' @param collapsed Start collapsed? (default TRUE)
#' @return HTML element for collapsible help
collapsible_help <- function(id, title, content, collapsed = TRUE) {
  collapse_class <- if (collapsed) "collapse" else "collapse show"
  button_class <- if (collapsed) "collapsed" else ""

  tags$div(
    class = "help-panel mb-3",
    tags$button(
      class = paste("btn btn-link text-decoration-none p-0", button_class),
      type = "button",
      `data-bs-toggle` = "collapse",
      `data-bs-target` = paste0("#", id),
      `aria-expanded` = if (collapsed) "false" else "true",
      `aria-controls` = id,
      style = "color: #0072B2; font-size: 0.9em;",
      tags$span(class = "me-1", HTML("&#9432;")),
      title,
      tags$span(class = "ms-1", style = "font-size: 0.8em;", HTML("&#9660;"))
    ),
    tags$div(
      id = id,
      class = collapse_class,
      style = "margin-top: 10px;",
      content
    )
  )
}

# =============================================================================
# Terminology Panel (Glossary)
# =============================================================================

#' Full terminology panel for help modal
terminology_panel_ui <- function() {
  tags$div(
    class = "terminology-panel",

    tags$h5("Terminology Guide", style = "margin-bottom: 15px;"),

    lapply(names(TERMINOLOGY), function(key) {
      term <- TERMINOLOGY[[key]]
      tags$div(
        style = "margin-bottom: 15px; padding-bottom: 15px; border-bottom: 1px solid #eee;",
        tags$h6(term$term, style = "color: #0072B2; margin-bottom: 5px;"),
        tags$p(term$full, style = "margin-bottom: 0; color: #555; font-size: 0.9em;")
      )
    })
  )
}

# =============================================================================
# JavaScript for Tooltips
# =============================================================================

#' JavaScript to initialize Bootstrap tooltips
tooltip_js <- function() {
  tags$script(HTML("
    $(document).ready(function() {
      // Initialize tooltips
      var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle=\"tooltip\"]'));
      var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
        return new bootstrap.Tooltip(tooltipTriggerEl);
      });

      // Reinitialize on Shiny content update
      $(document).on('shiny:value', function(event) {
        setTimeout(function() {
          var newTooltips = [].slice.call(document.querySelectorAll('[data-bs-toggle=\"tooltip\"]'));
          newTooltips.forEach(function(el) {
            if (!bootstrap.Tooltip.getInstance(el)) {
              new bootstrap.Tooltip(el);
            }
          });
        }, 100);
      });
    });
  "))
}
