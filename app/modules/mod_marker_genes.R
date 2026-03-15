# mod_marker_genes.R — Marker Genes module

marker_genes_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      title = "Options",
      width = 300,
      sliderInput(
        ns("n_markers"), "Top N markers per cell type:",
        min = 3, max = 20, value = 5, step = 1
      ),
      hr(),
      h6("Custom Gene Set"),
      selectizeInput(
        ns("custom_genes"), "Select genes for Dot Plot:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Type gene names...",
          maxOptions = 50,
          maxItems = 30
        )
      ),
      actionButton(ns("load_defaults"), "Load Default Markers",
                   class = "btn-sm btn-outline-primary"),
      hr(),
      downloadButton(ns("dl_heatmap"), "Download Heatmap (PDF)", class = "btn-sm download-btn"),
      downloadButton(ns("dl_dotplot"), "Download Dot Plot (PDF)", class = "btn-sm download-btn")
    ),
    layout_column_wrap(
      width = 1,
      card(
        full_screen = TRUE,
        card_header("Top Marker Genes Heatmap"),
        card_body(
          class = "plot-container",
          plotOutput(ns("heatmap"), height = "500px")
        )
      ),
      card(
        full_screen = TRUE,
        card_header("Dot Plot — Custom Gene Set"),
        card_body(
          class = "plot-container",
          plotOutput(ns("dotplot"), height = "500px")
        )
      )
    )
  )
}

marker_genes_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Server-side selectize
    updateSelectizeInput(session, "custom_genes", choices = gene_list,
                         server = TRUE)

    # Load default markers
    observeEvent(input$load_defaults, {
      defaults <- c("Myh4", "Myh1", "Myh2", "Myh7",
                     "Pax7", "Pecam1", "Pdgfra", "Ptprc",
                     "Rgs5", "Adipoq", "Chrna1")
      available <- defaults[defaults %in% gene_list]
      updateSelectizeInput(session, "custom_genes", selected = available)
    })

    # Top markers per cell type
    top_markers_df <- reactive({
      marker_genes %>%
        group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = input$n_markers) %>%
        ungroup()
    })

    # Heatmap
    heatmap_gg <- reactive({
      top_df <- top_markers_df()
      genes <- unique(top_df$gene)
      genes <- genes[genes %in% rownames(expr_matrix)]

      if (length(genes) == 0) return(NULL)

      # Compute mean expression per cell type
      ct_means <- do.call(rbind, lapply(cell_types, function(ct) {
        idx <- which(umap_meta$cell_type == ct)
        if (length(idx) == 0) return(rep(0, length(genes)))
        Matrix::rowMeans(expr_matrix[genes, idx, drop = FALSE])
      }))
      rownames(ct_means) <- cell_types
      colnames(ct_means) <- genes

      # Scale per gene
      ct_scaled <- scale(ct_means)
      ct_scaled[is.nan(ct_scaled)] <- 0

      # Convert to long format
      heat_df <- as.data.frame(ct_scaled) %>%
        mutate(cell_type = rownames(ct_scaled)) %>%
        pivot_longer(-cell_type, names_to = "gene", values_to = "scaled_expr")

      # Order genes by cluster
      gene_order <- top_df %>%
        arrange(cluster) %>%
        pull(gene) %>%
        unique()
      gene_order <- gene_order[gene_order %in% genes]

      heat_df$gene <- factor(heat_df$gene, levels = gene_order)
      heat_df$cell_type <- factor(heat_df$cell_type, levels = cell_types)

      ggplot(heat_df, aes(x = gene, y = cell_type, fill = scaled_expr)) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_viridis(option = "magma", name = "Scaled\nExpression") +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 10),
          panel.grid = element_blank()
        ) +
        labs(x = NULL, y = NULL)
    })

    output$heatmap <- renderPlot({
      p <- heatmap_gg()
      if (!is.null(p)) p
    })

    # Dot plot for custom genes
    dotplot_gg <- reactive({
      req(input$custom_genes)
      genes <- input$custom_genes[input$custom_genes %in% rownames(expr_matrix)]
      req(length(genes) > 0)

      # Compute per cell type stats
      dot_data <- do.call(rbind, lapply(cell_types, function(ct) {
        idx <- which(umap_meta$cell_type == ct)
        if (length(idx) == 0) return(NULL)

        do.call(rbind, lapply(genes, function(g) {
          vals <- expr_matrix[g, idx]
          data.frame(
            cell_type = ct,
            gene = g,
            avg_exp = mean(vals),
            pct_exp = mean(vals > 0) * 100,
            stringsAsFactors = FALSE
          )
        }))
      }))

      dot_data$cell_type <- factor(dot_data$cell_type, levels = rev(cell_types))
      dot_data$gene <- factor(dot_data$gene, levels = genes)

      ggplot(dot_data, aes(x = gene, y = cell_type)) +
        geom_point(aes(size = pct_exp, color = avg_exp)) +
        scale_color_viridis(option = "magma", name = "Avg\nExpression") +
        scale_size_continuous(name = "% Expressing", range = c(1, 10)) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_line(color = "grey90")
        ) +
        labs(x = NULL, y = NULL)
    })

    output$dotplot <- renderPlot({
      dotplot_gg()
    })

    # Downloads
    output$dl_heatmap <- downloadHandler(
      filename = function() "marker_heatmap.pdf",
      content = function(file) {
        p <- heatmap_gg()
        if (!is.null(p)) ggsave(file, p, width = 14, height = 8, device = "pdf")
      }
    )

    output$dl_dotplot <- downloadHandler(
      filename = function() "marker_dotplot.pdf",
      content = function(file) {
        ggsave(file, dotplot_gg(), width = 10, height = 8, device = "pdf")
      }
    )
  })
}
