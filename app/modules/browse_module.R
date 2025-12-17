# browse_module.R
# Shiny module for browsing and selecting tRNAs

# =============================================================================
# UI
# =============================================================================

browseUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Filters row
    fluidRow(
      column(4,
        selectInput(
          ns("amino_acid"),
          "Amino Acid",
          choices = NULL,  # Populated by server
          width = "100%"
        )
      ),
      column(4,
        selectInput(
          ns("compartment"),
          "Compartment",
          choices = c("All" = "all", "Nuclear" = "nuclear", "Mitochondrial" = "mitochondrial"),
          selected = "all",
          width = "100%"
        )
      ),
      column(4,
        tags$div(
          style = "padding-top: 25px;",
          actionButton(ns("select_all"), "Select All", class = "btn-sm"),
          actionButton(ns("clear_all"), "Clear All", class = "btn-sm btn-secondary")
        )
      )
    ),

    # Selection actions
    fluidRow(
      column(12,
        tags$div(
          style = "margin: 10px 0; padding: 10px; background-color: #f8f9fa; border-radius: 4px;",
          tags$span("Selected: "),
          tags$span(id = ns("selected_count"), "0"),
          tags$span(" tRNAs | "),
          actionButton(ns("mark_desired"), "Mark as Desired",
                       class = "btn-sm btn-desired"),
          actionButton(ns("mark_avoid"), "Mark as Avoid",
                       class = "btn-sm btn-avoid"),
          actionButton(ns("clear_selection"), "Clear Selection",
                       class = "btn-sm btn-outline-secondary")
        )
      )
    ),

    # Terminology guide (collapsible)
    fluidRow(
      column(12,
        tags$details(
          style = "margin-bottom: 15px;",
          tags$summary(
            style = "cursor: pointer; color: #0072B2;",
            "More about tRNA terminology"
          ),
          create_terminology_html()
        )
      )
    ),

    # Sequence browser
    fluidRow(
      column(12,
        uiOutput(ns("sequence_browser"))
      )
    ),

    # Hidden input to receive checkbox changes from JavaScript
    tags$input(
      type = "hidden",
      id = ns("checkbox_state"),
      value = ""
    ),

    # JavaScript for checkbox handling (including group checkboxes)
    tags$script(HTML(sprintf("
      // Handle group checkbox - select/deselect all members
      $(document).on('change', '#%s .group-checkbox', function() {
        var members = $(this).data('group-members');
        if (members) {
          var memberIds = members.split(',');
          var isChecked = $(this).is(':checked');
          // Check/uncheck all member checkboxes
          memberIds.forEach(function(id) {
            $('input.member-checkbox[data-id=\"' + id + '\"]').prop('checked', isChecked);
          });
        }
        updateCheckedState();
      });

      // Handle member checkbox
      $(document).on('change', '#%s .member-checkbox', function() {
        updateCheckedState();
      });

      // Handle regular tRNA row checkbox
      $(document).on('change', '#%s .trna-row .trna-checkbox', function() {
        updateCheckedState();
      });

      // Collect all checked IDs and send to Shiny
      function updateCheckedState() {
        var checked = [];

        // Get checked from regular rows
        $('#%s .trna-row .trna-checkbox:checked').each(function() {
          var row = $(this).closest('.trna-row');
          var id = row.data('id');
          if (id) checked.push(id);
        });

        // Get checked from group members
        $('#%s .member-checkbox:checked').each(function() {
          var id = $(this).data('id');
          if (id && checked.indexOf(id) === -1) checked.push(id);
        });

        // Get members of checked groups
        $('#%s .group-checkbox:checked').each(function() {
          var members = $(this).data('group-members');
          if (members) {
            members.split(',').forEach(function(id) {
              if (checked.indexOf(id) === -1) checked.push(id);
            });
          }
        });

        Shiny.setInputValue('%s', checked.join(','), {priority: 'event'});
      }

      // Toggle expand icon on collapse
      $(document).on('show.bs.collapse', '#%s .group-members', function() {
        $(this).prev().find('.expand-icon').html('\\u25BC');
      });
      $(document).on('hide.bs.collapse', '#%s .group-members', function() {
        $(this).prev().find('.expand-icon').html('\\u25B6');
      });
    ", ns("sequence_browser"), ns("sequence_browser"), ns("sequence_browser"),
       ns("sequence_browser"), ns("sequence_browser"), ns("sequence_browser"),
       ns("checkbox_state"), ns("sequence_browser"), ns("sequence_browser"))))
  )
}

# =============================================================================
# Server
# =============================================================================

browseServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Local reactive values for this module
    local <- reactiveValues(
      checked_ids = character()  # Currently checked (not yet assigned to desired/avoid)
    )

    # -------------------------------------------------------------------------
    # Update amino acid choices when data loads
    # -------------------------------------------------------------------------
    observe({
      req(values$trna_data)

      # Get unique amino acids
      aa_list <- sort(unique(values$trna_data$amino_acid))

      # Format choices with full names
      choices <- setNames(aa_list, sapply(aa_list, format_amino_acid))

      updateSelectInput(session, "amino_acid", choices = choices)
    })

    # -------------------------------------------------------------------------
    # Filter data based on selections
    # -------------------------------------------------------------------------
    filtered_data <- reactive({
      req(values$trna_data, input$amino_acid)

      df <- values$trna_data

      # Filter by amino acid
      df <- df[df$amino_acid == input$amino_acid, ]

      # Filter by compartment
      if (input$compartment != "all") {
        df <- df[df$compartment == input$compartment | is.na(df$compartment), ]
      }

      df
    })

    # -------------------------------------------------------------------------
    # Render sequence browser
    # -------------------------------------------------------------------------
    output$sequence_browser <- renderUI({
      req(filtered_data())

      df <- filtered_data()

      if (nrow(df) == 0) {
        return(tags$p("No tRNAs found for this selection."))
      }

      # Get organism name
      org_name <- names(ORGANISMS)[ORGANISMS == isolate(input$organism)]
      if (length(org_name) == 0) org_name <- NULL

      # Render the view
      render_amino_acid_view(
        df,
        amino_acid = input$amino_acid,
        selected_ids = values$selected_ids,
        selection_types = values$selection_types,
        organism_name = org_name
      )
    })

    # -------------------------------------------------------------------------
    # Handle checkbox changes
    # -------------------------------------------------------------------------
    observeEvent(input$checkbox_state, {
      if (is.null(input$checkbox_state) || input$checkbox_state == "") {
        local$checked_ids <- character()
      } else {
        local$checked_ids <- strsplit(input$checkbox_state, ",")[[1]]
      }

      # Update count display
      session$sendCustomMessage("updateCount", list(
        id = ns("selected_count"),
        count = length(local$checked_ids)
      ))
    }, ignoreNULL = FALSE)

    # JavaScript to update count
    tags$script(HTML("
      Shiny.addCustomMessageHandler('updateCount', function(data) {
        var el = document.getElementById(data.id);
        if (el) el.innerText = data.count;
      });
    "))

    # -------------------------------------------------------------------------
    # Selection actions
    # -------------------------------------------------------------------------

    # Select all visible
    observeEvent(input$select_all, {
      df <- filtered_data()
      if (!is.null(df) && nrow(df) > 0) {
        # Use JavaScript to check all boxes
        session$sendCustomMessage("checkAll", ns("sequence_browser"))
      }
    })

    # Clear all visible
    observeEvent(input$clear_all, {
      session$sendCustomMessage("uncheckAll", ns("sequence_browser"))
    })

    # JavaScript handlers for check/uncheck all
    tags$script(HTML("
      Shiny.addCustomMessageHandler('checkAll', function(containerId) {
        $('#' + containerId + ' .trna-checkbox').prop('checked', true).trigger('change');
      });
      Shiny.addCustomMessageHandler('uncheckAll', function(containerId) {
        $('#' + containerId + ' .trna-checkbox').prop('checked', false).trigger('change');
      });
    "))

    # Mark as desired
    observeEvent(input$mark_desired, {
      req(length(local$checked_ids) > 0)

      # Add to selection types
      new_types <- values$selection_types
      for (id in local$checked_ids) {
        new_types[[id]] <- "desired"
      }
      values$selection_types <- new_types

      # Update selected_ids
      values$selected_ids <- unique(c(values$selected_ids, local$checked_ids))

      showNotification(
        paste0("Marked ", length(local$checked_ids), " tRNAs as desired targets"),
        type = "message",
        duration = 2
      )
    })

    # Mark as avoid
    observeEvent(input$mark_avoid, {
      req(length(local$checked_ids) > 0)

      # Add to selection types
      new_types <- values$selection_types
      for (id in local$checked_ids) {
        new_types[[id]] <- "avoid"
      }
      values$selection_types <- new_types

      # Update selected_ids
      values$selected_ids <- unique(c(values$selected_ids, local$checked_ids))

      showNotification(
        paste0("Marked ", length(local$checked_ids), " tRNAs to avoid"),
        type = "message",
        duration = 2
      )
    })

    # Clear selection (remove from desired/avoid)
    observeEvent(input$clear_selection, {
      req(length(local$checked_ids) > 0)

      # Remove from selection types
      new_types <- values$selection_types
      for (id in local$checked_ids) {
        new_types[[id]] <- NULL
      }
      values$selection_types <- new_types

      # Remove from selected_ids
      values$selected_ids <- setdiff(values$selected_ids, local$checked_ids)

      showNotification(
        paste0("Cleared selection for ", length(local$checked_ids), " tRNAs"),
        type = "message",
        duration = 2
      )
    })

    # -------------------------------------------------------------------------
    # Return selection state
    # -------------------------------------------------------------------------
    return(list(
      selected_ids = reactive(values$selected_ids),
      selection_types = reactive(values$selection_types)
    ))

  })
}
