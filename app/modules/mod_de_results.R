# mod_de_results.R — Differential Expression Results module

de_results_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      title = "DE Filters",
      width = 300,
      selectInput(
        ns("cell_type"), "Cell Type:",
        choices = NULL
      ),
      sliderInput(
        ns("logfc_thresh"), "|log2FC| threshold:",
        min = 0, max = 3, value = 0.25, step = 0.05
      ),
      sliderInput(
        ns("pval_thresh"), "-log10(padj) threshold:",
        min = 0, max = 50, value = 1.3, step = 0.1
      ),
      hr(),
      downloadButton(ns("dl_table"), "Download Table (CSV)", class = "btn-sm download-btn"),
      downloadButton(ns("dl_volcano"), "Download Volcano (PDF)", class = "btn-sm download-btn")
    ),
    layout_column_wrap(
      width = 1,
      card(
        card_header("Volcano Plot"),
        card_body(
          class = "plot-container",
          plotlyOutput(ns("volcano"), height = "450px")
        )
      ),
      card(
        card_header("DE Genes Table"),
        card_body(
          DTOutput(ns("de_table"))
        )
      )
    )
  )
}

de_results_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Populate cell type choices
    observe({
      ct_choices <- names(de_results)
      updateSelectInput(session, "cell_type", choices = ct_choices,
                        selected = ct_choices[1])
    })

    # Reactive: filtered DE data
    de_data <- reactive({
      req(input$cell_type)
      req(input$cell_type %in% names(de_results))

      df <- de_results[[input$cell_type]]
      df$neg_log10_padj <- -log10(pmax(df$p_val_adj, 1e-300))

      # Classify significance
      df$significance <- "NS"
      df$significance[abs(df$avg_log2FC) > input$logfc_thresh &
                        df$neg_log10_padj > input$pval_thresh] <- "Significant"
      df$significance[abs(df$avg_log2FC) > input$logfc_thresh &
                        df$neg_log10_padj > input$pval_thresh &
                        df$avg_log2FC > 0] <- "Up in KIC"
      df$significance[abs(df$avg_log2FC) > input$logfc_thresh &
                        df$neg_log10_padj > input$pval_thresh &
                        df$avg_log2FC < 0] <- "Down in KIC"
      df
    })

    # Volcano plot (plotly)
    volcano_gg <- reactive({
      df <- de_data()

      sig_colors <- c("Up in KIC" = "#E74C3C", "Down in KIC" = "#3498DB", "NS" = "#BDC3C7")

      p <- ggplot(df, aes(x = avg_log2FC, y = neg_log10_padj,
                           color = significance, text = gene)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = sig_colors, name = NULL) +
        geom_vline(xintercept = c(-input$logfc_thresh, input$logfc_thresh),
                   linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = input$pval_thresh,
                   linetype = "dashed", color = "grey50") +
        theme_minimal(base_size = 14) +
        labs(
          x = "log2 Fold Change (KIC vs Control)",
          y = "-log10(adjusted p-value)",
          title = paste(input$cell_type, "— KIC vs Control")
        ) +
        theme(legend.position = "top")
      p
    })

    output$volcano <- renderPlotly({
      p <- volcano_gg()
      ggplotly(p, tooltip = c("text", "x", "y")) %>%
        layout(legend = list(orientation = "h", y = 1.1))
    })

    # DE table
    output$de_table <- renderDT({
      df <- de_data()
      sig_df <- df %>%
        filter(significance != "NS") %>%
        select(gene, avg_log2FC, pct.1, pct.2, p_val_adj, significance) %>%
        arrange(p_val_adj) %>%
        mutate(
          avg_log2FC = round(avg_log2FC, 3),
          pct.1 = round(pct.1, 3),
          pct.2 = round(pct.2, 3),
          p_val_adj = signif(p_val_adj, 3)
        )

      datatable(
        sig_df,
        colnames = c("Gene", "log2FC", "% KIC", "% Control", "Adj. P-value", "Direction"),
        options = list(pageLength = 15, scrollX = TRUE),
        selection = "single",
        rownames = FALSE
      )
    })

    # Downloads
    output$dl_table <- downloadHandler(
      filename = function() {
        ct_clean <- gsub("[/ ]", "_", input$cell_type)
        paste0("DE_", ct_clean, ".csv")
      },
      content = function(file) {
        df <- de_data() %>%
          filter(significance != "NS") %>%
          select(gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj, significance) %>%
          arrange(p_val_adj)
        write.csv(df, file, row.names = FALSE)
      }
    )

    output$dl_volcano <- downloadHandler(
      filename = function() {
        ct_clean <- gsub("[/ ]", "_", input$cell_type)
        paste0("volcano_", ct_clean, ".pdf")
      },
      content = function(file) {
        ggsave(file, volcano_gg(), width = 8, height = 6, device = "pdf")
      }
    )
  })
}
