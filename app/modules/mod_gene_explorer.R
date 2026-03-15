# mod_gene_explorer.R — Gene Explorer module

gene_explorer_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      title = "Gene Search",
      width = 300,
      div(class = "gene-search",
        selectizeInput(
          ns("gene"), "Select Gene:",
          choices = NULL,
          options = list(
            placeholder = "Type gene name...",
            maxOptions = 50
          )
        )
      ),
      radioButtons(
        ns("plot_type"), "Expression Plot:",
        choices = c("Violin" = "violin", "Dot Plot" = "dotplot"),
        selected = "violin"
      ),
      checkboxInput(ns("split_umap"), "Split UMAP by Condition", value = TRUE),
      hr(),
      downloadButton(ns("dl_umap"), "Download UMAP", class = "btn-sm download-btn"),
      downloadButton(ns("dl_expr"), "Download Expression Plot", class = "btn-sm download-btn")
    ),
    layout_column_wrap(
      width = 1,
      card(
        card_header(textOutput(ns("umap_title"))),
        card_body(
          class = "plot-container",
          plotOutput(ns("umap_plot"), height = "450px")
        )
      ),
      card(
        card_header(textOutput(ns("expr_title"))),
        card_body(
          class = "plot-container",
          plotOutput(ns("expr_plot"), height = "400px")
        )
      )
    )
  )
}

gene_explorer_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Server-side gene selectize (case-insensitive)
    updateSelectizeInput(session, "gene", choices = gene_list,
                         server = TRUE,
                         selected = character(0))

    # Reactive: get expression for selected gene
    gene_expr <- reactive({
      req(input$gene)
      req(input$gene %in% rownames(expr_matrix))

      vals <- expr_matrix[input$gene, ]
      data.frame(
        UMAP1 = umap_meta$UMAP1,
        UMAP2 = umap_meta$UMAP2,
        expression = as.numeric(vals),
        cell_type = umap_meta$cell_type,
        condition = umap_meta$condition,
        stringsAsFactors = FALSE
      )
    })

    # Titles
    output$umap_title <- renderText({
      if (is.null(input$gene) || input$gene == "") {
        "UMAP — Select a gene"
      } else {
        paste0("UMAP — ", input$gene, " Expression")
      }
    })

    output$expr_title <- renderText({
      if (is.null(input$gene) || input$gene == "") {
        "Expression Plot — Select a gene"
      } else {
        paste0(input$gene, " — ", ifelse(input$plot_type == "violin", "Violin", "Dot"), " Plot")
      }
    })

    # UMAP plot
    umap_gg <- reactive({
      df <- gene_expr()
      # Sort so high expression on top
      df <- df[order(df$expression), ]

      p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = expression)) +
        geom_point(size = 0.3, alpha = 0.7) +
        scale_color_viridis(option = "magma", name = input$gene) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "right",
          panel.grid = element_blank(),
          axis.title = element_text(size = 12)
        ) +
        labs(x = "UMAP 1", y = "UMAP 2")

      if (input$split_umap) {
        p <- p + facet_wrap(~condition)
      }
      p
    })

    output$umap_plot <- renderPlot({
      req(input$gene)
      umap_gg()
    })

    # Expression plot (violin or dotplot)
    expr_gg <- reactive({
      df <- gene_expr()

      if (input$plot_type == "violin") {
        ggplot(df, aes(x = cell_type, y = expression, fill = condition)) +
          geom_violin(scale = "width", trim = TRUE, alpha = 0.8,
                      position = position_dodge(width = 0.8)) +
          scale_fill_manual(values = c("Control" = "#3498DB", "KIC" = "#E74C3C")) +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            legend.position = "top"
          ) +
          labs(x = NULL, y = "Normalized Expression", fill = "Condition")
      } else {
        # Dot plot: percent expressing + mean expression per cell type + condition
        dot_df <- df %>%
          group_by(cell_type, condition) %>%
          summarise(
            avg_exp = mean(expression),
            pct_exp = mean(expression > 0) * 100,
            .groups = "drop"
          )

        ggplot(dot_df, aes(x = condition, y = cell_type)) +
          geom_point(aes(size = pct_exp, color = avg_exp)) +
          scale_color_viridis(option = "magma", name = "Avg Expression") +
          scale_size_continuous(name = "% Expressing", range = c(1, 10)) +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.y = element_text(size = 11),
            legend.position = "right"
          ) +
          labs(x = NULL, y = NULL)
      }
    })

    output$expr_plot <- renderPlot({
      req(input$gene)
      expr_gg()
    })

    # Downloads
    output$dl_umap <- downloadHandler(
      filename = function() paste0(input$gene, "_umap.pdf"),
      content = function(file) {
        ggsave(file, umap_gg(), width = 10, height = 5, device = "pdf")
      }
    )

    output$dl_expr <- downloadHandler(
      filename = function() paste0(input$gene, "_expression.pdf"),
      content = function(file) {
        ggsave(file, expr_gg(), width = 10, height = 6, device = "pdf")
      }
    )
  })
}
