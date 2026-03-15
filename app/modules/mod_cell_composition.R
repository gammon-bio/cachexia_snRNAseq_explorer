# mod_cell_composition.R — Cell Composition module

cell_composition_ui <- function(id) {
  ns <- NS(id)

  layout_column_wrap(
    width = 1 / 2,
    heights_equal = "row",
    card(
      card_header("UMAP by Cell Type"),
      card_body(
        class = "plot-container",
        plotlyOutput(ns("umap_celltype"), height = "500px")
      )
    ),
    card(
      card_header("UMAP by Condition"),
      card_body(
        class = "plot-container",
        plotlyOutput(ns("umap_condition"), height = "500px")
      )
    ),
    card(
      full_screen = TRUE,
      card_header("Cell Type Proportions"),
      card_body(
        class = "plot-container",
        plotOutput(ns("bar_chart"), height = "400px")
      ),
      card_footer(
        downloadButton(ns("dl_bar"), "Download (PDF)", class = "btn-sm")
      )
    ),
    card(
      card_header("Cell Counts"),
      card_body(
        DTOutput(ns("count_table"))
      )
    )
  )
}

cell_composition_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Subsample for plotly performance
    umap_sub <- reactive({
      n <- nrow(umap_meta)
      if (n > 20000) {
        idx <- sample(n, 20000)
        umap_meta[idx, ]
      } else {
        umap_meta
      }
    })

    # UMAP by cell type (plotly)
    output$umap_celltype <- renderPlotly({
      df <- umap_sub()
      p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = cell_type,
                           text = cell_type)) +
        geom_point(size = 0.3, alpha = 0.7) +
        scale_color_manual(values = color_palette, name = "Cell Type") +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(), legend.position = "right") +
        labs(x = "UMAP 1", y = "UMAP 2")

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(itemsizing = "constant"))
    })

    # UMAP by condition (plotly)
    output$umap_condition <- renderPlotly({
      df <- umap_sub()
      p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = condition,
                           text = condition)) +
        geom_point(size = 0.3, alpha = 0.7) +
        scale_color_manual(values = c("Control" = "#3498DB", "KIC" = "#E74C3C")) +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(), legend.position = "right") +
        labs(x = "UMAP 1", y = "UMAP 2")

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(itemsizing = "constant"))
    })

    # Stacked bar chart
    bar_gg <- reactive({
      df <- umap_meta %>%
        count(condition, cell_type) %>%
        group_by(condition) %>%
        mutate(proportion = n / sum(n)) %>%
        ungroup()

      ggplot(df, aes(x = condition, y = proportion, fill = cell_type)) +
        geom_col(position = "fill", width = 0.6) +
        scale_fill_manual(values = color_palette, name = "Cell Type") +
        scale_y_continuous(labels = percent_format()) +
        theme_minimal(base_size = 14) +
        labs(x = NULL, y = "Proportion", title = "Cell Type Proportions by Condition") +
        theme(legend.position = "right")
    })

    output$bar_chart <- renderPlot({
      bar_gg()
    })

    # Count table
    output$count_table <- renderDT({
      df <- summary_stats %>%
        select(cell_type, condition, n_cells) %>%
        pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0) %>%
        mutate(Total = Control + KIC) %>%
        arrange(desc(Total)) %>%
        rename(`Cell Type` = cell_type)

      datatable(df, options = list(pageLength = 15, dom = "t"),
                rownames = FALSE)
    })

    # Download
    output$dl_bar <- downloadHandler(
      filename = function() "cell_composition.pdf",
      content = function(file) {
        ggsave(file, bar_gg(), width = 8, height = 6, device = "pdf")
      }
    )
  })
}
