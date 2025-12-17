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
# Sequence Grouping and Comparison
# =============================================================================

#' Group tRNAs by sequence identity
#'
#' Groups tRNAs with identical sequences together and identifies differences
#' between non-identical sequences within the same isoacceptor.
#'
#' @param trna_df Data frame with tRNA data (must have 'sequence' and 'id' columns)
#' @return Data frame with added columns: sequence_group, group_size, is_representative,
#'         diff_from_ref, diff_positions
#' @export
group_by_sequence <- function(trna_df) {
  if (nrow(trna_df) == 0) return(trna_df)

  # Assign sequence groups (identical sequences get same group)
  unique_seqs <- unique(trna_df$sequence)
  seq_to_group <- setNames(seq_along(unique_seqs), unique_seqs)

  trna_df$sequence_group <- seq_to_group[trna_df$sequence]

  # Count group sizes
  group_counts <- table(trna_df$sequence_group)
  trna_df$group_size <- as.integer(group_counts[as.character(trna_df$sequence_group)])

  # Mark first in each group as representative
  trna_df$is_representative <- !duplicated(trna_df$sequence_group)

  # Use first sequence as reference for comparisons
  ref_seq <- trna_df$sequence[1]

  # Calculate differences from reference
  trna_df$diff_from_ref <- mapply(count_differences, trna_df$sequence,
                                   MoreArgs = list(ref = ref_seq))
  trna_df$diff_positions <- mapply(find_difference_positions, trna_df$sequence,
                                    MoreArgs = list(ref = ref_seq), SIMPLIFY = FALSE)

  trna_df
}

#' Count nucleotide differences between two sequences
#'
#' @param seq1 First sequence
#' @param seq2 Second sequence (reference)
#' @return Integer count of differences
count_differences <- function(seq1, ref) {
  if (nchar(seq1) != nchar(ref)) {
    return(abs(nchar(seq1) - nchar(ref)) +
           sum(strsplit(seq1, "")[[1]][1:min(nchar(seq1), nchar(ref))] !=
               strsplit(ref, "")[[1]][1:min(nchar(seq1), nchar(ref))]))
  }
  sum(strsplit(seq1, "")[[1]] != strsplit(ref, "")[[1]])
}

#' Find positions where two sequences differ
#'
#' @param seq1 First sequence
#' @param ref Reference sequence
#' @return Character vector describing differences (e.g., "pos 45: G→A")
find_difference_positions <- function(seq1, ref) {
  if (seq1 == ref) return(character(0))

  chars1 <- strsplit(seq1, "")[[1]]
  chars_ref <- strsplit(ref, "")[[1]]

  # Handle length differences
  min_len <- min(length(chars1), length(chars_ref))

  diffs <- character(0)

  # Check each position
  for (i in seq_len(min_len)) {
    if (chars1[i] != chars_ref[i]) {
      diffs <- c(diffs, paste0(i, ":", chars_ref[i], "\u2192", chars1[i]))
    }
  }

  # Note length differences
  if (length(chars1) > length(chars_ref)) {
    diffs <- c(diffs, paste0("+", length(chars1) - length(chars_ref), "nt at 3'"))
  } else if (length(chars1) < length(chars_ref)) {
    diffs <- c(diffs, paste0("-", length(chars_ref) - length(chars1), "nt at 3'"))
  }

  diffs
}

#' Format difference summary for display
#'
#' Shows % identity for many differences, or specific positions for few (1-5)
#'
#' @param n_diff Number of differences
#' @param diff_positions Character vector from find_difference_positions
#' @param seq_length Length of the sequence
#' @return Formatted string
format_differences <- function(n_diff, diff_positions = NULL, seq_length = 75) {
  if (n_diff == 0) return("")

  # Calculate % identity
  pct_identity <- round(100 * (seq_length - n_diff) / seq_length, 0)

  if (n_diff <= 5 && !is.null(diff_positions) && length(diff_positions) > 0) {
    # Show specific positions for small number of diffs
    paste0(pct_identity, "% identity (", paste(diff_positions, collapse = ", "), ")")
  } else {
    # Just show % identity for many diffs
    paste0(pct_identity, "% identity (", n_diff, " nt differ)")
  }
}

# =============================================================================
# Multiple Sequence Alignment View
# =============================================================================

#' Find variable positions across a set of sequences
#'
#' @param sequences Character vector of sequences
#' @return Integer vector of 1-indexed positions where sequences differ
find_variable_positions <- function(sequences) {
  if (length(sequences) <= 1) return(integer(0))


  # Get max length

max_len <- max(nchar(sequences))

  # Pad sequences to same length
  padded <- sapply(sequences, function(s) {
    if (nchar(s) < max_len) {
      paste0(s, paste(rep("-", max_len - nchar(s)), collapse = ""))
    } else {
      s
    }
  })

  # Split into character matrices
  char_mat <- do.call(rbind, strsplit(padded, ""))

  # Find positions where not all characters are the same
  variable <- sapply(seq_len(ncol(char_mat)), function(i) {
    length(unique(char_mat[, i])) > 1
  })

  which(variable)
}

#' Compute consensus sequence
#'
#' @param sequences Character vector of sequences
#' @return Character string with most common base at each position
compute_consensus <- function(sequences) {
  max_len <- max(nchar(sequences))

  # Pad sequences
 padded <- sapply(sequences, function(s) {
    if (nchar(s) < max_len) {
      paste0(s, paste(rep("-", max_len - nchar(s)), collapse = ""))
    } else {
      s
    }
  })

  # Get consensus (most common) at each position
  char_mat <- do.call(rbind, strsplit(padded, ""))

  sapply(seq_len(ncol(char_mat)), function(i) {
    tab <- table(char_mat[, i])
    names(tab)[which.max(tab)]
  }) |> paste(collapse = "")
}

#' Collapse identical sequences into groups
#'
#' @param group_df Data frame with tRNA data
#' @return List with unique_seqs (data frame) and group_info (list of member IDs per unique seq)
collapse_identical_sequences <- function(group_df) {
  # Group by sequence
  seq_groups <- split(group_df, group_df$sequence)

  # For each unique sequence, keep first row as representative
  unique_df <- do.call(rbind, lapply(seq_groups, function(g) g[1, ]))

  # Store member IDs for each group
  group_info <- lapply(seq_groups, function(g) g$id)
  names(group_info) <- unique_df$sequence

  list(
    unique_df = unique_df,
    group_info = group_info
  )
}

#' Get color for nucleotide variant (Sanger-style with Okabe-Ito palette)
#'
#' @param base Character, the nucleotide
#' @return CSS color string
get_base_color <- function(base) {
  switch(toupper(base),
    "A" = "#009E73",  # Green (Okabe-Ito bluish green)
    "T" = "#D55E00",  # Red (Okabe-Ito vermillion)
    "U" = "#D55E00",  # Red (same as T)
    "G" = "#F0E442",  # Yellow (Okabe-Ito)
    "C" = "#0072B2",  # Blue (Okabe-Ito blue)
    "#999999"         # Gray for gaps/unknown
  )
}

#' Render MSA view with consensus and collapsed identical sequences
#'
#' @param group_df Data frame with tRNA data for one isoacceptor
#' @param selected_ids Character vector of selected tRNA IDs
#' @param selection_types Named list mapping IDs to "desired" or "avoid"
#' @return HTML element with MSA visualization
render_msa_view <- function(group_df, selected_ids, selection_types) {
  if (nrow(group_df) == 0) return(NULL)

  # Collapse identical sequences
  collapsed <- collapse_identical_sequences(group_df)
  unique_df <- collapsed$unique_df
  group_info <- collapsed$group_info

  sequences <- unique_df$sequence
  n_unique <- length(sequences)

  # Compute consensus
  consensus <- compute_consensus(sequences)
  consensus_chars <- strsplit(consensus, "")[[1]]
  max_len <- nchar(consensus)

  # Get anticodon positions
  ac_start <- unique_df$anticodon_start[1]
  ac_end <- unique_df$anticodon_end[1]

  # Build the MSA display
  tags$div(
    class = "msa-container",
    style = "overflow-x: auto; font-family: monospace; font-size: 0.85em;",

    # Consensus row
    tags$div(
      class = "msa-row consensus-row",
      style = "display: flex; align-items: center; padding: 4px; margin-bottom: 4px; background-color: #f0f0f0; border-bottom: 2px solid #ccc;",

      # Spacer for checkbox alignment (matches checkbox + margin)
      tags$span(style = "width: 21px; flex-shrink: 0;"),

      tags$span(
        style = "width: 190px; flex-shrink: 0; font-weight: bold; color: #333; font-size: 0.9em;",
        "Consensus"
      ),

      {
        seq_html <- sapply(seq_len(max_len), function(pos) {
          char <- consensus_chars[pos]
          is_anticodon <- !is.na(ac_start) && pos >= ac_start && pos <= ac_end

          style <- if (is_anticodon) {
            paste0("background-color: ", COLORS$anticodon, "; font-weight: bold;")
          } else {
            "color: #333; font-weight: bold;"
          }

          sprintf('<span style="%s">%s</span>', style, char)
        })

        tags$span(
          class = "msa-sequence",
          style = "white-space: pre;",
          HTML(paste(seq_html, collapse = ""))
        )
      }
    ),

    # Each unique sequence row
    lapply(seq_len(n_unique), function(i) {
      row <- unique_df[i, ]
      seq_chars <- strsplit(row$sequence, "")[[1]]
      seq_len_i <- length(seq_chars)

      # Get all member IDs for this sequence
      member_ids <- group_info[[row$sequence]]
      n_members <- length(member_ids)

      # Check selection status of members
      any_selected <- any(member_ids %in% selected_ids)
      member_sel_types <- sapply(member_ids, function(id) selection_types[[id]])
      n_desired <- sum(member_sel_types == "desired", na.rm = TRUE)
      n_avoid <- sum(member_sel_types == "avoid", na.rm = TRUE)

      # Selection styling
      sel_style <- if (n_desired > 0 && n_avoid == 0) {
        paste0("border-left: 3px solid ", COLORS$desired, ";")
      } else if (n_avoid > 0 && n_desired == 0) {
        paste0("border-left: 3px solid ", COLORS$avoid, ";")
      } else if (n_desired > 0 && n_avoid > 0) {
        "border-left: 3px solid #999;"  # Mixed
      } else ""

      # Unique ID for collapsible (if multi-gene)
      collapse_id <- paste0("genes-", gsub("[^a-zA-Z0-9]", "-", row$id))

      # Container for row + collapsible
      tags$div(
        class = "sequence-group-container",

        # Main row
        tags$div(
          class = "msa-row",
          style = paste0(
            "display: flex; align-items: center; padding: 2px 4px; margin: 1px 0; ",
            sel_style,
            if (any_selected) "background-color: #fafafa;" else ""
          ),
          `data-id` = row$id,
          `data-members` = paste(member_ids, collapse = ","),

          # Checkbox
          tags$input(
            type = "checkbox",
            class = "trna-checkbox group-checkbox",
            `data-group-members` = paste(member_ids, collapse = ","),
            checked = if (any_selected) "checked" else NULL,
            style = "margin-right: 8px; flex-shrink: 0;"
          ),

          # Label with count and expand button for multi-gene groups
          tags$span(
            class = "msa-id",
            style = "width: 190px; flex-shrink: 0; font-size: 0.9em;",
            if (n_members > 1) {
              tags$span(
                style = "cursor: pointer; color: #666;",
                `data-bs-toggle` = "collapse",
                `data-bs-target` = paste0("#", collapse_id),
                tags$span(class = "expand-icon", style = "margin-right: 4px;", "\u25B6"),
                paste0(n_members, " identical genes")
              )
            } else {
              tags$span(style = "color: #666;", member_ids[1])
            }
          ),

          # Sequence showing only differences from consensus
          {
            seq_html <- sapply(seq_len(max_len), function(pos) {
              char <- if (pos <= seq_len_i) seq_chars[pos] else "-"
              cons_char <- consensus_chars[pos]
              is_anticodon <- !is.na(ac_start) && pos >= ac_start && pos <= ac_end
              is_diff <- (char != cons_char)

              if (is_anticodon) {
                style <- paste0("background-color: ", COLORS$anticodon, "; font-weight: bold;")
                sprintf('<span style="%s">%s</span>', style, char)
              } else if (is_diff) {
                color <- get_base_color(char)
                # Use dark text for yellow (G), white for others
                text_color <- if (toupper(char) == "G") "#333" else "white"
                style <- paste0("background-color: ", color, "; color: ", text_color, "; font-weight: bold;")
                sprintf('<span style="%s">%s</span>', style, char)
              } else {
                sprintf('<span style="color: #999;">%s</span>', char)
              }
            })

            tags$span(
              class = "msa-sequence",
              style = "white-space: pre;",
              HTML(paste(seq_html, collapse = ""))
            )
          }
        ),

        # Collapsible member list for multi-gene groups
        if (n_members > 1) {
          tags$div(
            id = collapse_id,
            class = "collapse",
            style = "margin-left: 30px; padding: 5px 10px; background-color: #f9f9f9; border-left: 2px solid #ddd; font-size: 0.85em;",
            tags$div(
              style = "color: #555; margin-bottom: 5px; font-weight: 500;",
              "Member genes:"
            ),
            lapply(member_ids, function(mid) {
              sel_type <- selection_types[[mid]]
              tags$div(
                style = "padding: 2px 0;",
                tags$input(
                  type = "checkbox",
                  class = "trna-checkbox member-checkbox",
                  `data-id` = mid,
                  checked = if (mid %in% selected_ids) "checked" else NULL,
                  style = "margin-right: 6px;"
                ),
                tags$span(style = "font-family: monospace;", mid),
                if (!is.null(sel_type)) {
                  if (sel_type == "desired") {
                    tags$span(style = paste0("color: ", COLORS$desired, "; margin-left: 5px;"), "\u2713")
                  } else {
                    tags$span(style = paste0("color: ", COLORS$avoid, "; margin-left: 5px;"), "\u2717")
                  }
                }
              )
            })
          )
        }
      )
    }),

    # Legend
    tags$div(
      class = "msa-legend",
      style = "margin-top: 10px; font-size: 0.8em; color: #666; display: flex; flex-wrap: wrap; gap: 10px; align-items: center;",
      tags$span("Variants: "),
      tags$span(style = "background-color: #009E73; color: white; padding: 1px 6px;", "A"),
      tags$span(style = "background-color: #D55E00; color: white; padding: 1px 6px;", "T"),
      tags$span(style = "background-color: #F0E442; color: #333; padding: 1px 6px;", "G"),
      tags$span(style = "background-color: #0072B2; color: white; padding: 1px 6px;", "C"),
      tags$span(style = "margin-left: 10px;", "|"),
      tags$span(style = paste0("background-color: ", COLORS$anticodon, "; padding: 1px 6px;"), "Anticodon")
    )
  )
}

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

#' Render rows grouped by sequence identity
#'
#' Shows identical sequences as collapsed groups with expandable member lists
#' and displays differences for non-identical sequences.
#'
#' @param group_df Data frame with sequence grouping (from group_by_sequence)
#' @param selected_ids Character vector of selected tRNA IDs
#' @param selection_types Named list mapping IDs to "desired" or "avoid"
#' @return List of HTML elements
render_grouped_rows <- function(group_df, selected_ids, selection_types) {
  # Process each unique sequence group
  unique_groups <- unique(group_df$sequence_group)
  n_unique <- length(unique_groups)

  lapply(seq_along(unique_groups), function(idx) {
    grp_id <- unique_groups[idx]
    members <- group_df[group_df$sequence_group == grp_id, ]
    n_members <- nrow(members)
    representative <- members[1, ]

    # Check if any members are selected
    member_ids <- members$id
    any_selected <- any(member_ids %in% selected_ids)

    # Get selection types for members
    member_sel_types <- lapply(member_ids, function(id) selection_types[[id]])

    # Determine label: first group vs isodecoder (different sequence)
    is_first <- (idx == 1)

    if (n_members == 1) {
      # Single sequence
      seq_label <- if (!is_first && n_unique > 1) {
        tags$span(
          class = "isodecoder-label",
          style = "margin-left: 10px; font-size: 0.8em; color: #E69F00; font-style: italic;",
          "isodecoder"
        )
      } else {
        NULL
      }

      tags$div(
        class = "sequence-group single",
        format_trna_row_html(
          representative,
          selected = representative$id %in% selected_ids,
          selection_type = selection_types[[representative$id]],
          checkbox_id = paste0("chk-", gsub("[^a-zA-Z0-9]", "-", representative$id))
        ),
        seq_label
      )

    } else {
      # Multiple identical sequences - show as collapsible group
      group_id <- paste0("grp-", gsub("[^a-zA-Z0-9]", "-", representative$id))

      # Determine group selection state
      n_desired <- sum(unlist(member_sel_types) == "desired", na.rm = TRUE)
      n_avoid <- sum(unlist(member_sel_types) == "avoid", na.rm = TRUE)

      group_status <- if (n_desired == n_members) {
        tags$span(style = paste0("color: ", COLORS$desired, ";"), " \u2713 all desired")
      } else if (n_avoid == n_members) {
        tags$span(style = paste0("color: ", COLORS$avoid, ";"), " \u2717 all avoid")
      } else if (n_desired > 0 || n_avoid > 0) {
        tags$span(style = "color: #666;",
                  paste0(" (", n_desired, " desired, ", n_avoid, " avoid)"))
      } else {
        NULL
      }

      # Label for isodecoder groups (different sequence from first)
      seq_label <- if (!is_first && n_unique > 1) {
        tags$span(
          class = "isodecoder-label",
          style = "font-size: 0.8em; color: #E69F00; font-style: italic; margin-left: 5px;",
          "isodecoder"
        )
      } else {
        NULL
      }

      tags$div(
        class = "sequence-group multiple",
        style = "margin: 4px 0;",

        # Collapsed view with expand button
        tags$div(
          class = "group-header",
          style = paste0(
            "display: flex; align-items: center; padding: 4px 8px; ",
            "background-color: ", if (any_selected) "#f5f5f5" else "#fafafa", "; ",
            "border-radius: 4px; cursor: pointer; ",
            "white-space: nowrap; overflow-x: auto;"
          ),
          `data-bs-toggle` = "collapse",
          `data-bs-target` = paste0("#", group_id),

          # Checkbox for group (selects all)
          tags$input(
            type = "checkbox",
            class = "trna-checkbox group-checkbox",
            id = paste0("chk-", group_id),
            `data-group-members` = paste(member_ids, collapse = ","),
            checked = if (any_selected) "checked" else NULL,
            style = "margin-right: 10px;"
          ),

          # Expand icon
          tags$span(
            class = "expand-icon",
            style = "margin-right: 8px; color: #666;",
            "\u25B6"  # Right-pointing triangle
          ),

          # Sequence with anticodon
          format_sequence_html(
            representative$sequence,
            representative$anticodon_start,
            representative$anticodon_end
          ),

          # Group count badge
          tags$span(
            class = "badge bg-secondary",
            style = "margin-left: 10px;",
            paste0(n_members, " identical")
          ),

          seq_label,
          group_status
        ),

        # Expandable member list
        tags$div(
          id = group_id,
          class = "collapse group-members",
          style = "padding-left: 30px; border-left: 2px solid #e0e0e0; margin-left: 15px;",

          lapply(seq_len(n_members), function(i) {
            member <- members[i, ]
            tags$div(
              class = "group-member",
              style = "padding: 2px 0; font-size: 0.9em;",
              tags$input(
                type = "checkbox",
                class = "trna-checkbox member-checkbox",
                `data-id` = member$id,
                checked = if (member$id %in% selected_ids) "checked" else NULL,
                style = "margin-right: 8px;"
              ),
              tags$span(
                class = "trna-id",
                style = "color: #666;",
                member$id
              ),
              if (!is.null(selection_types[[member$id]])) {
                if (selection_types[[member$id]] == "desired") {
                  tags$span(style = paste0("color: ", COLORS$desired, "; margin-left: 5px;"), "\u2713")
                } else {
                  tags$span(style = paste0("color: ", COLORS$avoid, "; margin-left: 5px;"), "\u2717")
                }
              }
            )
          })
        )
      )
    }
  })
}

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

      # Count unique sequences
      n_unique <- length(unique(group$sequence))

      tags$div(
        class = "isoacceptor-group",
        style = "margin-bottom: 20px;",

        # Group header with unique count
        tags$div(
          class = "isoacceptor-header",
          style = paste0(
            "font-weight: bold; ",
            "padding: 5px 10px; ",
            "background-color: #f0f0f0; ",
            "border-left: 4px solid ", COLORS$header, ";"
          ),
          paste0("Isoacceptor: ", amino_acid, "-", anticodon,
                 " (", n_genes, " gene", if (n_genes != 1) "s" else "", ")"),
          if (n_unique < n_genes) {
            tags$span(
              style = "font-weight: normal; color: #666; margin-left: 10px;",
              paste0("\u2014 ", n_unique, " unique sequence",
                     if (n_unique != 1) "s" else "")
            )
          }
        ),

        # MSA view showing aligned sequences with differences highlighted
        tags$div(
          class = "isoacceptor-msa",
          style = "padding: 10px 0;",
          render_msa_view(group, selected_ids, selection_types)
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

    .sequence-group {
      margin: 2px 0;
    }

    .group-header {
      white-space: nowrap;
      overflow-x: auto;
    }

    .group-header:hover {
      background-color: #f0f0f0 !important;
    }

    .diff-info {
      flex-shrink: 0;
    }

    .isoacceptor-rows {
      overflow-x: auto;
    }
  ")))
}
