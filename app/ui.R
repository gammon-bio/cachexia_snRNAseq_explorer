# ui.R
# bslib page_navbar layout with 5 tabs

ui <- page_navbar(
  title = "Zhang et al. Cachexia snRNAseq Explorer",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    "navbar-bg" = "#2C3E50"
  ),
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    useWaiter()
  ),

  # Tab 1: Gene Explorer
  nav_panel(
    title = "Gene Explorer",
    icon = icon("dna"),
    gene_explorer_ui("gene_explorer")
  ),

  # Tab 2: Differential Expression
  nav_panel(
    title = "Differential Expression",
    icon = icon("chart-bar"),
    de_results_ui("de_results")
  ),

  # Tab 3: Cell Composition
  nav_panel(
    title = "Cell Composition",
    icon = icon("circle-nodes"),
    cell_composition_ui("cell_composition")
  ),

  # Tab 4: Marker Genes
  nav_panel(
    title = "Marker Genes",
    icon = icon("list"),
    marker_genes_ui("marker_genes")
  ),

  # Tab 5: About
  nav_panel(
    title = "About",
    icon = icon("info-circle"),
    layout_column_wrap(
      width = 1,
      fill = FALSE,
      card(
        card_header("About This App"),
        card_body(
          h4("Zhang et al. Cachexia snRNAseq Explorer"),
          p("This interactive application allows exploration of single-nucleus
             RNA-seq data from mouse skeletal muscle comparing ",
            tags$b("control"), " vs ", tags$b("KIC cachexia"), " model mice."),
          hr(),
          h5("Paper"),
          p("Zhang et al. ",
            tags$em("A myogenin-myostatin pathway drives cancer cachexia-induced
                     muscle atrophy."),
            " Cell Reports, 2024."),
          p(tags$b("PMID: "), "39116208"),
          p(tags$b("GEO: "),
            tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272083",
                   target = "_blank", "GSE272083"),
            " (snRNA-seq: GSE272085)"),
          hr(),
          h5("Methods"),
          tags$ul(
            tags$li("Single-nucleus multiome (snRNA-seq + snATAC-seq) on mouse
                     gastrocnemius-plantaris skeletal muscle"),
            tags$li("This app shows snRNA-seq (gene expression) data only"),
            tags$li("Normalization: SCTransform v2 (used for clustering, UMAP,
                     and expression visualization in this app)"),
            tags$li("Clustering: Seurat v5 (Louvain algorithm, resolution 0.5, 30 PCs)"),
            tags$li("DE testing: Wilcoxon rank-sum test on log-normalized RNA assay
                     (standard Seurat NormalizeData; used instead of SCT for
                     cross-condition DE to avoid multi-model correction issues)"),
            tags$li("9,379 nuclei total (3,914 control + 5,465 KIC)")
          ),
          hr(),
          h5("Dataset Summary"),
          tableOutput("about_summary_table"),
          hr(),
          h5("Important Note"),
          p(tags$em("This dataset represents pooled samples without biological
                     replicates at the library level. DE results should be
                     interpreted as exploratory and validated experimentally."),
            class = "text-muted")
        )
      )
    )
  )
)
