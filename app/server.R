# server.R

server <- function(input, output, session) {

  # Initialize waiter
  w <- Waiter$new(
    html = spin_fading_circles(),
    color = "rgba(44, 62, 80, 0.7)"
  )

  # Gene Explorer module
  gene_explorer_server("gene_explorer")

  # DE Results module
  de_results_server("de_results")

  # Cell Composition module
  cell_composition_server("cell_composition")

  # Marker Genes module
  marker_genes_server("marker_genes")

  # About tab: summary table
  output$about_summary_table <- renderTable({
    summary_stats %>%
      group_by(cell_type) %>%
      summarise(
        Control = sum(n_cells[condition == "Control"]),
        KIC = sum(n_cells[condition == "KIC"]),
        Total = sum(n_cells),
        .groups = "drop"
      ) %>%
      arrange(desc(Total)) %>%
      rename(`Cell Type` = cell_type)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
}
