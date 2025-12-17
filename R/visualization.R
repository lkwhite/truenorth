# visualization.R
# HTML visualization functions for tRNA sequence display
# Includes anticodon detection, sequence formatting, and educational content

library(htmltools)
library(dplyr)

# =============================================================================
# Color Palette - Okabe-Ito colorblind-friendly colors
# =============================================================================

# Okabe-Ito color palette (colorblind-friendly)
OKABE_ITO <- list(
  orange = "#E69F00",
  sky_blue = "#56B4E9",
  bluish_green = "#009E73",
  yellow = "#F0E442",
  blue = "#0072B2",
  vermillion = "#D55E00",
  reddish_purple = "#CC79A7",
  black = "#000000"
)

# Semantic color assignments
COLORS <- list(
  anticodon = OKABE_ITO$sky_blue,      # Highlight anticodon
  desired = OKABE_ITO$bluish_green,    # Desired targets (green-ish)
  avoid = OKABE_ITO$vermillion,        # Avoid targets (orange-red)
  header = OKABE_ITO$blue,             # Section headers
  info = OKABE_ITO$sky_blue,           # Info boxes
  warning = OKABE_ITO$orange           # Warning/tip boxes
)

# =============================================================================
# Anticodon Detection
# =============================================================================

#' Find anticodon position in a tRNA sequence
#'
#' Locates the anticodon by searching for the expected 3-mer from the header
#' in the expected position range. Handles U→T conversion.
#'
#' @param sequence Character string of tRNA sequence
#' @param anticodon_from_header Anticodon from tRNA header (may use U, e.g., "UGG")
#' @param min_pos Minimum expected position (default 30)
#' @param max_pos Maximum expected position (default 45)
#' @return List with start, end, anticodon, or NULL if not found
#' @export
find_anticodon_position <- function(sequence, anticodon_from_header,
                                    min_pos = 30, max_pos = 45) {
  if (is.na(sequence) || is.na(anticodon_from_header)) {
    return(NULL)
  }

  # Convert U to T in the anticodon (header may use RNA notation)
  anticodon_dna <- gsub("U", "T", toupper(anticodon_from_header))

  # Find all occurrences of the anticodon in the sequence
  seq_upper <- toupper(sequence)
  matches <- gregexpr(anticodon_dna, seq_upper, fixed = TRUE)[[1]]

  # No matches found

  if (matches[1] == -1) {
    return(NULL)
  }

  # Filter to expected position range
  valid_matches <- matches[matches >= min_pos & matches <= max_pos]

  if (length(valid_matches) == 0) {
    # Check if there are matches outside the expected range
    if (length(matches) > 0) {
      # Return the closest match to position 35, but flag it
      closest <- matches[which.min(abs(matches - 35))]
      return(list(
        start = closest,
        end = closest + 2,
        anticodon = anticodon_dna,
        validated = FALSE,
        note = paste0("Anticodon found at position ", closest,
                      " (outside expected range ", min_pos, "-", max_pos, ")")
      ))
    }
    return(NULL)
  }

  # If multiple matches in range, prefer one closest to position 35
  best_match <- valid_matches[which.min(abs(valid_matches - 35))]

  list(
    start = best_match,
    end = best_match + 2,
    anticodon = anticodon_dna,
    validated = TRUE,
    note = NULL
  )
}

#' Add anticodon position to tRNA data frame
#'
#' Computes anticodon positions for all tRNAs in a data frame
#'
#' @param trna_df Data frame from load_organism_data()
#' @return Data frame with anticodon_start, anticodon_end columns added
#' @export
add_anticodon_positions <- function(trna_df) {
  # Apply find_anticodon_position to each row
  positions <- mapply(
    function(seq, ac) {
      pos <- find_anticodon_position(seq, ac)
      if (is.null(pos)) {
        c(NA_integer_, NA_integer_, FALSE)
      } else {
        c(pos$start, pos$end, pos$validated)
      }
    },
    trna_df$sequence,
    trna_df$anticodon,
    SIMPLIFY = TRUE
  )

  trna_df$anticodon_start <- as.integer(positions[1, ])
  trna_df$anticodon_end <- as.integer(positions[2, ])
  trna_df$anticodon_validated <- as.logical(positions[3, ])

  trna_df
}

# =============================================================================
# HTML Sequence Formatting
# =============================================================================

#' Format a tRNA sequence as HTML with anticodon highlighted
#'
#' @param sequence Character string of tRNA sequence
#' @param anticodon_start Start position of anticodon (1-indexed)
#' @param anticodon_end End position of anticodon (1-indexed, inclusive)
#' @param highlight_color CSS color for anticodon (default: Okabe-Ito sky blue)
#' @return HTML string with highlighted anticodon
#' @export
format_sequence_html <- function(sequence, anticodon_start, anticodon_end,
                                 highlight_color = COLORS$anticodon) {
  if (is.na(anticodon_start) || is.na(anticodon_end)) {
    # No anticodon position - return plain sequence
    return(tags$span(class = "trna-sequence", sequence))
  }

  # Split sequence into parts
  before <- substr(sequence, 1, anticodon_start - 1)
  anticodon <- substr(sequence, anticodon_start, anticodon_end)
  after <- substr(sequence, anticodon_end + 1, nchar(sequence))

  # Build HTML with highlighted anticodon
  tags$span(
    class = "trna-sequence",
    style = "font-family: monospace;",
    before,
    tags$span(
      class = "anticodon",
      style = paste0("background-color: ", highlight_color, "; ",
                     "font-weight: bold; ",
                     "padding: 1px 2px; ",
                     "border-radius: 2px;"),
      anticodon
    ),
    after
  )
}

#' Format a single tRNA row for display
#'
#' Creates a complete row with checkbox, sequence, and ID
#'
#' @param trna_row Single row from tRNA data frame
#' @param selected Logical, is this tRNA selected?
#' @param selection_type "desired", "avoid", or NULL
#' @param checkbox_id Optional ID for checkbox element
#' @return HTML div element
#' @export
format_trna_row_html <- function(trna_row, selected = FALSE,
                                 selection_type = NULL, checkbox_id = NULL) {
  # Selection indicator styling (using Okabe-Ito colors)
  selection_style <- ""
  selection_icon <- ""
  if (!is.null(selection_type)) {
    if (selection_type == "desired") {
      selection_style <- paste0("border-left: 3px solid ", COLORS$desired, ";")
      selection_icon <- tags$span(
        style = paste0("color: ", COLORS$desired, "; margin-right: 5px;"),
        "\u2713"  # checkmark
      )
    } else if (selection_type == "avoid") {
      selection_style <- paste0("border-left: 3px solid ", COLORS$avoid, ";")
      selection_icon <- tags$span(
        style = paste0("color: ", COLORS$avoid, "; margin-right: 5px;"),
        "\u2717"  # X mark
      )
    }
  }

  # Build the row
  tags$div(
    class = paste("trna-row", if (selected) "selected" else ""),
    style = paste0(
      "display: flex; ",
      "align-items: center; ",
      "padding: 4px 8px; ",
      "margin: 2px 0; ",
      "background-color: ", if (selected) "#f0f0f0" else "transparent", "; ",
      selection_style
    ),
    `data-id` = trna_row$id,

    # Checkbox
    tags$input(
      type = "checkbox",
      class = "trna-checkbox",
      id = checkbox_id,
      checked = if (selected) "checked" else NULL,
      style = "margin-right: 10px;"
    ),

    # Selection icon
    selection_icon,

    # Sequence with highlighted anticodon
    format_sequence_html(
      trna_row$sequence,
      trna_row$anticodon_start,
      trna_row$anticodon_end
    ),

    # tRNA ID
    tags$span(
      class = "trna-id",
      style = "margin-left: 15px; color: #666; font-size: 0.9em;",
      trna_row$id
    )
  )
}

# =============================================================================
# Isoacceptor Family View
# =============================================================================

#' Render all tRNAs for an amino acid, grouped by isoacceptor
#'
#' Creates the main visualization showing all tRNAs for an amino acid,
#' organized by anticodon (isoacceptor), with educational context.
#'
#' @param trna_df Data frame from load_organism_data() (with anticodon positions)
#' @param amino_acid Three-letter amino acid code to display
#' @param selected_ids Character vector of selected tRNA IDs
#' @param selection_types Named list mapping IDs to "desired" or "avoid"
#' @param organism_name Display name for organism (e.g., "Human")
#' @return HTML div element containing the full view
#' @export
render_amino_acid_view <- function(trna_df, amino_acid,
                                   selected_ids = character(0),
                                   selection_types = list(),
                                   organism_name = NULL) {
  # Filter to amino acid
  aa_trnas <- trna_df[trna_df$amino_acid == amino_acid, ]

  if (nrow(aa_trnas) == 0) {
    return(tags$div(
      class = "amino-acid-view empty",
      tags$p("No tRNAs found for ", amino_acid)
    ))
  }

  # Get full amino acid name
  aa_names <- c(
    Ala = "Alanine", Arg = "Arginine", Asn = "Asparagine", Asp = "Aspartic acid",
    Cys = "Cysteine", Gln = "Glutamine", Glu = "Glutamic acid", Gly = "Glycine",
    His = "Histidine", Ile = "Isoleucine", Leu = "Leucine", Lys = "Lysine",
    Met = "Methionine", Phe = "Phenylalanine", Pro = "Proline", Ser = "Serine",
    Thr = "Threonine", Trp = "Tryptophan", Tyr = "Tyrosine", Val = "Valine",
    iMet = "Initiator Methionine", fMet = "Formyl-Methionine",
    SeC = "Selenocysteine", Sup = "Suppressor"
  )
  aa_full <- aa_names[amino_acid]
  if (is.na(aa_full)) aa_full <- amino_acid

  # Group by isoacceptor (anticodon)
  isoacceptors <- split(aa_trnas, aa_trnas$anticodon)

  # Build the view
  tags$div(
    class = "amino-acid-view",
    style = "font-family: sans-serif;",

    # Header
    tags$div(
      class = "aa-header",
      style = paste0(
        "padding: 10px; ",
        "background-color: #f8f9fa; ",
        "border-bottom: 2px solid #dee2e6; ",
        "margin-bottom: 15px;"
      ),
      tags$h3(
        style = "margin: 0;",
        paste0(aa_full, " (", amino_acid, ") tRNAs"),
        if (!is.null(organism_name)) paste0(" - ", organism_name),
        tags$span(
          style = "float: right; font-size: 0.8em; color: #666;",
          paste0("[", nrow(aa_trnas), " total]")
        )
      )
    ),

    # Terminology mini-guide
    tags$div(
      class = "terminology-guide",
      style = paste0(
        "background-color: #e7f3ff; ",
        "border: 1px solid #b3d7ff; ",
        "border-radius: 4px; ",
        "padding: 10px; ",
        "margin-bottom: 15px; ",
        "font-size: 0.9em;"
      ),
      tags$strong("Terminology: "),
      tags$span(
        "Isoacceptor = same anticodon (e.g., all ", amino_acid, "-",
        names(isoacceptors)[1], "). ",
        "Gene family = copies from same ancestral gene (", amino_acid, "-",
        names(isoacceptors)[1], "-1-* are all family 1)."
      )
    ),

    # Isoacceptor groups
    lapply(names(isoacceptors), function(anticodon) {
      group <- isoacceptors[[anticodon]]
      n_genes <- nrow(group)

      # Sort by gene family and copy number
      group <- group[order(group$gene_family, group$copy_number), ]

      tags$div(
        class = "isoacceptor-group",
        style = "margin-bottom: 20px;",

        # Group header
        tags$div(
          class = "isoacceptor-header",
          style = paste0(
            "font-weight: bold; ",
            "padding: 5px 10px; ",
            "background-color: #f0f0f0; ",
            "border-left: 4px solid ", COLORS$header, ";"
          ),
          paste0("Isoacceptor: ", amino_acid, "-", anticodon,
                 " (", n_genes, " gene", if (n_genes != 1) "s" else "", ")")
        ),

        # tRNA rows
        tags$div(
          class = "isoacceptor-rows",
          style = "padding-left: 10px;",
          lapply(seq_len(nrow(group)), function(i) {
            row <- group[i, ]
            is_selected <- row$id %in% selected_ids
            sel_type <- selection_types[[row$id]]
            format_trna_row_html(
              row,
              selected = is_selected,
              selection_type = sel_type,
              checkbox_id = paste0("chk-", gsub("[^a-zA-Z0-9]", "-", row$id))
            )
          })
        )
      )
    }),

    # Selection summary footer
    tags$div(
      class = "selection-summary",
      style = paste0(
        "margin-top: 15px; ",
        "padding: 10px; ",
        "background-color: #f8f9fa; ",
        "border-top: 1px solid #dee2e6;"
      ),
      {
        n_desired <- if (length(selection_types) == 0) 0 else sum(unlist(selection_types) == "desired")
        n_avoid <- if (length(selection_types) == 0) 0 else sum(unlist(selection_types) == "avoid")
        tags$span(
          paste0("Selection: ", n_desired, " desired, ", n_avoid, " avoid")
        )
      }
    )
  )
}

# =============================================================================
# Educational Content
# =============================================================================

#' Create HTML for the terminology education panel
#'
#' Returns a standalone HTML panel explaining tRNA nomenclature
#'
#' @return HTML div element
#' @export
create_terminology_html <- function() {
  tags$div(
    class = "terminology-panel",
    style = paste0(
      "font-family: sans-serif; ",
      "padding: 15px; ",
      "background-color: #f8f9fa; ",
      "border: 1px solid #dee2e6; ",
      "border-radius: 4px; ",
      "line-height: 1.6;"
    ),

    tags$h4(
      style = "margin-top: 0;",
      "Understanding tRNA Nomenclature"
    ),

    # Example breakdown
    tags$div(
      style = "font-family: monospace; margin: 15px 0; background: white; padding: 10px; border-radius: 4px;",
      tags$pre(
        style = "margin: 0; font-size: 0.95em;",
        "Example: tRNA-Ala-AGC-12-3\n",
        "              |   |   |  |\n",
        "              |   |   |  +-- Copy number (3rd copy)\n",
        "              |   |   +----- Gene family (family 12)\n",
        "              |   +--------- Anticodon (AGC)\n",
        "              +------------- Amino acid (Alanine)"
      )
    ),

    # Definitions
    tags$dl(
      style = "margin: 0;",

      tags$dt(style = "font-weight: bold; margin-top: 10px;", "ISOACCEPTOR"),
      tags$dd(
        style = "margin-left: 20px;",
        "All tRNAs with the same anticodon (e.g., all Ala-AGC tRNAs). ",
        "These decode the same codon(s)."
      ),

      tags$dt(style = "font-weight: bold; margin-top: 10px;", "ISODECODER"),
      tags$dd(
        style = "margin-left: 20px;",
        "tRNAs with the same anticodon but different body sequences. ",
        tags$span(style = "color: #dc3545;", "Important for probe specificity!"),
        " Two isodecoders may not cross-hybridize to the same probe."
      ),

      tags$dt(style = "font-weight: bold; margin-top: 10px;", "GENE FAMILY"),
      tags$dd(
        style = "margin-left: 20px;",
        "Copies from the same ancestral gene (e.g., Ala-AGC-1-1, Ala-AGC-1-2, Ala-AGC-1-3 ",
        "are all family 1). Usually have identical or nearly identical sequences."
      )
    ),

    # Practical tip
    tags$div(
      style = "margin-top: 15px; padding: 10px; background: #fff3cd; border-radius: 4px; font-size: 0.9em;",
      tags$strong("Probe Design Tip: "),
      "To hit an entire gene family, design a probe against any member. ",
      "To discriminate between isodecoders, look for regions where their sequences differ."
    )
  )
}

# =============================================================================
# Utility Functions
# =============================================================================

#' Get CSS styles for tRNA visualization
#'
#' Returns CSS stylesheet for use in Shiny apps (using Okabe-Ito colors)
#'
#' @return HTML style tag
#' @export
get_visualization_css <- function() {
  tags$style(HTML(paste0("
    .trna-sequence {
      font-family: 'Courier New', monospace;
      font-size: 0.95em;
      letter-spacing: 1px;
    }

    .anticodon {
      background-color: ", COLORS$anticodon, ";
      font-weight: bold;
      padding: 1px 2px;
      border-radius: 2px;
    }

    .trna-row {
      cursor: pointer;
      transition: background-color 0.15s;
    }

    .trna-row:hover {
      background-color: #e9ecef !important;
    }

    .trna-row.selected {
      background-color: #f0f0f0;
    }

    .isoacceptor-group {
      margin-bottom: 20px;
    }

    .isoacceptor-header {
      font-weight: bold;
      padding: 5px 10px;
      background-color: #f0f0f0;
      border-left: 4px solid ", COLORS$header, ";
    }

    .amino-acid-view {
      max-width: 1200px;
    }

    .terminology-guide {
      background-color: ", paste0(COLORS$info, "22"), ";
      border: 1px solid ", COLORS$info, ";
      border-radius: 4px;
      padding: 10px;
      margin-bottom: 15px;
      font-size: 0.9em;
    }
  ")))
}
