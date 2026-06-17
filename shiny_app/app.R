library(shiny)
library(bslib)
library(bsicons)
library(dplyr)
library(readr)
library(plotly)
library(DT)

# --- Safe data loading ---
safe_read <- function(path) {
  tryCatch(
    read_csv(path, show_col_types = FALSE),
    error = function(e) { message("Failed to load: ", path); NULL }
  )
}

mr_screen    <- safe_read("data/STable1_full_MR_screen.csv")
master_ev    <- safe_read("data/STable_master_evidence.csv")
druggability <- safe_read("data/STable5_druggability.csv")
mediation    <- safe_read("data/STable4_mediation_integrated_evidence.csv")

# Format volcano data
if (!is.null(mr_screen)) {
  mr_screen <- mr_screen %>%
    mutate(
      log_or      = log(OR),
      neg_log_p   = -log10(pvalue),
      hit_label   = if_else(FDR < 0.05, "FDR < 0.05", "Non-significant"),
      hover_text  = paste0(
        "<b>", protein, "</b>",
        "<br>Cancer: ", gsub("_GCST.*", "", cancer_outcome),
        "<br>OR: ",  round(OR, 3),
        "<br>P: ",   signif(pvalue, 3),
        "<br>FDR: ", signif(FDR, 3)
      )
    )
}

cancer_choices <- if (!is.null(mr_screen)) {
  v  <- unique(mr_screen$cancer_outcome)
  nm <- gsub("_GCST.*", "", v)
  setNames(v, nm)
} else character(0)

# Tier colour palette (used in two tables)
tier_pal <- c(T1 = "#dff0d8", T2a = "#d9edf7", T2b = "#fcf8e3", T2c = "#f5f5f5")

# --- Theme ---
app_theme <- bs_theme(
  version   = 5,
  bootswatch = "zephyr",
  primary   = "#005b96",
  secondary = "#b3cde0",
  success   = "#28a745",
  info      = "#17a2b8"
)

# Reusable static-image helper
static_img <- function(src, ...) {
  div(style = "text-align:center;",
    img(src = src, style = paste(
      "max-width:100%; height:auto; border-radius:8px;",
      "box-shadow:0 4px 8px rgba(0,0,0,.12);", ...
    ))
  )
}

# --- UI ---
ui <- page_navbar(
  title   = "Multi-Omic Network MR",
  theme   = app_theme,
  bg      = "#005b96",
  inverse = TRUE,

  # 1. Overview ----------------------------------------------------------------
  nav_panel(
    title = "Overview", icon = bs_icon("info-circle"),
    page_sidebar(
      sidebar = sidebar(
        title = "About",
        p("Interactive supplement to:"),
        strong("Multi-omic triangulation of circulating proteins identifies novel breast cancer causal candidates."),
        hr(),
        em("Modhukur V et al. 2026"),
        hr(),
        tags$a(bs_icon("github"), " GitHub",
               href   = "https://github.com/vijayachitrabio/multiomic-network-mr-cancer",
               target = "_blank")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header(class = "bg-primary text-white", "Key findings"),
          card_body(
            tags$ul(
              tags$li("701 circulating proteins screened; cis-pQTL instruments from FinnGen Olink (N = 619)."),
              tags$li("17 protein–cancer associations survived BH-FDR correction."),
              tags$li(strong("5 Tier 1 (MR + dual coloc):"),
                      " EFNA1, TNFRSF6B, ATRAID, FGF5, ABO."),
              tags$li(strong("3 Tier 2a (MR + ABF-only coloc):"),
                      " SNX15, PM20D1, UMOD."),
              tags$li(strong("2 Tier 2b (MR + partial coloc):"),
                      " APOE, TSPAN8."),
              tags$li(strong("7 Tier 2c (MR-supported):"),
                      " remaining proteins.")
            )
          )
        ),
        card(
          card_header(class = "bg-secondary", "Study design"),
          card_body(static_img("sfig_mr_design_multiomic.png"))
        )
      )
    )
  ),

  # 2. MR Screen ---------------------------------------------------------------
  nav_panel(
    title = "MR Screen", icon = bs_icon("graph-up"),
    card(
      full_screen = TRUE,
      card_header(
        div(class = "d-flex justify-content-between align-items-center",
            "Proteome-wide MR screen — interactive volcano",
            selectInput("cancer_filter", "Cancer:",
                        choices = c("All cancers" = "all", cancer_choices),
                        width = "240px")
        )
      ),
      card_body(
        p(class = "text-muted small mb-1",
          tags$span(style = "color:#d9534f; font-weight:bold;", "●"),
          " Coloured = FDR < 0.05 (BH). X-axis = ln(OR) = β. ",
          "No Bonferroni line — significance defined by FDR across all tests."),
        plotlyOutput("volcano_plot", height = "570px")
      )
    )
  ),

  # 3. Forest Plot -------------------------------------------------------------
  nav_panel(
    title = "Forest Plot", icon = bs_icon("bar-chart-steps"),
    card(
      full_screen = TRUE,
      card_header("Forest plot — 17 FDR-significant protein–cancer associations"),
      card_body(plotlyOutput("forest_plot", height = "600px"))
    )
  ),

  # 4. Colocalization ----------------------------------------------------------
  nav_panel(
    title = "Colocalization", icon = bs_icon("diagram-2"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("coloc.abf vs coloc.SuSiE — method comparison"),
        card_body(
          div(style = "max-height:520px; overflow-y:auto;",
            static_img("fig_coloc_method_comparison.png"))
        )
      ),
      card(
        card_header("PPH4 results per protein"),
        card_body(DTOutput("coloc_table"))
      )
    )
  ),

  # 5. ER Subtypes -------------------------------------------------------------
  nav_panel(
    title = "ER Subtypes", icon = bs_icon("gender-female"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("ER-subtype forest plot (from manuscript)"),
        card_body(static_img("fig2_er_subtype_forest.png"))
      ),
      card(
        card_header("ER+ vs ER− OR concordance (interactive)"),
        card_body(plotlyOutput("er_scatter", height = "430px")),
        card_footer(class = "text-muted small",
          "Diagonal = perfect concordance. Dashed lines = OR = 1.")
      )
    )
  ),

  # 6. MAGMA -------------------------------------------------------------------
  nav_panel(
    title = "MAGMA", icon = bs_icon("bar-chart"),
    card(
      full_screen = TRUE,
      card_header("MAGMA gene-level associations for 17 MR-prioritised protein genes"),
      card_body(static_img("fig6_magma_mr_hits.png"))
    )
  ),

  # 7. Replication -------------------------------------------------------------
  nav_panel(
    title = "Replication", icon = bs_icon("check2-circle"),
    navset_card_tab(
      nav_panel("ARIC SomaScan",
        static_img("fig7_aric_replication.png")
      ),
      nav_panel("OpenGWAS",
        static_img("fig8_opengwas_replication.png")
      )
    )
  ),

  # 8. Immune Profiling --------------------------------------------------------
  nav_panel(
    title = "Immune", icon = bs_icon("shield-check"),
    navset_card_tab(
      nav_panel("TCGA-BRCA",
        static_img("fig9_tcga_immune_correlations.png")
      ),
      nav_panel("CPTAC-BRCA",
        static_img("fig10_cptac_protein_immune_correlations.png")
      ),
      nav_panel("TISCH scRNA-seq",
        div(style = "text-align:center; overflow-y:auto; max-height:720px; padding:12px;",
          h6("Major lineage"),
          static_img("fig11_tisch_scrna_majorlineage_heatmap.png",
                     "margin-bottom:20px;"),
          h6("Minor lineage"),
          static_img("fig12_tisch_scrna_minorlineage_heatmap.png")
        )
      )
    )
  ),

  # 9. Evidence Table ----------------------------------------------------------
  nav_panel(
    title = "Evidence", icon = bs_icon("table"),
    card(
      full_screen = TRUE,
      card_header("Multi-layer evidence — 17 FDR-significant candidates"),
      card_body(DTOutput("evidence_table"))
    )
  ),

  # 10. Druggability -----------------------------------------------------------
  nav_panel(
    title = "Druggability", icon = bs_icon("capsule"),
    card(
      full_screen = TRUE,
      card_header("Open Targets druggability — MR-prioritised proteins"),
      card_body(DTOutput("druggability_table"))
    )
  ),

  # 11. Mediation --------------------------------------------------------------
  nav_panel(
    title = "Mediation", icon = bs_icon("bezier2"),
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Mediation paths: protein → metabolite → cancer"),
        card_body(
          div(style = "max-height:520px; overflow-y:auto;",
            static_img("fig4_mediation_paths.png"))
        )
      ),
      card(
        card_header("Mediation evidence table"),
        card_body(DTOutput("mediation_table"))
      )
    )
  ),

  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode"))
)

# --- Server ------------------------------------------------------------------
server <- function(input, output, session) {

  # Reactive plotly layout that respects dark mode
  pl_layout <- reactive({
    if (isTRUE(input$dark_mode)) {
      list(paper_bgcolor = "#1a1a2e", plot_bgcolor = "#1a1a2e",
           font = list(color = "#e0e0e0"), gridcolor = "#333")
    } else {
      list(paper_bgcolor = "white", plot_bgcolor = "white",
           font = list(color = "#333"), gridcolor = "#eee")
    }
  })

  # Volcano plot
  output$volcano_plot <- renderPlotly({
    req(mr_screen)
    d <- mr_screen
    if (!is.null(input$cancer_filter) && input$cancer_filter != "all")
      d <- filter(d, cancer_outcome == input$cancer_filter)

    bl <- pl_layout()
    plot_ly(d,
            x = ~log_or, y = ~neg_log_p,
            type = "scatter", mode = "markers",
            color = ~hit_label,
            colors = c("FDR < 0.05" = "#d9534f", "Non-significant" = "#bbbbbb"),
            text = ~hover_text, hoverinfo = "text",
            marker = list(size = 8, opacity = 0.75,
                          line = list(width = 0.5, color = "white"))) %>%
      layout(
        xaxis  = list(title = "ln(OR) = β", zeroline = TRUE,
                      zerolinecolor = "#999", gridcolor = bl$gridcolor),
        yaxis  = list(title = "−log₁₀(p)", zeroline = FALSE,
                      gridcolor = bl$gridcolor),
        legend = list(title = list(text = "Significance")),
        paper_bgcolor = bl$paper_bgcolor,
        plot_bgcolor  = bl$plot_bgcolor,
        font          = bl$font
      )
  })

  # Forest plot
  output$forest_plot <- renderPlotly({
    req(master_ev)
    fp <- master_ev %>%
      filter(!is.na(mr_or)) %>%
      arrange(tier_short, cancer_mr, mr_or) %>%
      mutate(
        label = paste0(protein, " (", tier_short, " — ", cancer_mr, ")"),
        hover = paste0(
          "<b>", protein, "</b><br>",
          "Cancer: ", cancer_mr, "<br>",
          "Tier: ", tier_short, "<br>",
          "OR: ", round(mr_or, 3),
          " (", round(mr_or_lo, 3), "–", round(mr_or_hi, 3), ")<br>",
          "FDR: ", signif(mr_fdr, 3)
        )
      )
    fp$label <- factor(fp$label, levels = fp$label)

    bl <- pl_layout()
    plot_ly(fp, x = ~mr_or, y = ~label,
            type = "scatter", mode = "markers",
            color = ~cancer_mr,
            error_x = list(type = "data", symmetric = FALSE,
                           array = ~mr_or_hi - mr_or,
                           arrayminus = ~mr_or - mr_or_lo),
            text = ~hover, hoverinfo = "text",
            marker = list(size = 10)) %>%
      layout(
        xaxis  = list(title = "Odds ratio (95% CI)", gridcolor = bl$gridcolor),
        yaxis  = list(title = "", autorange = "reversed"),
        shapes = list(list(type = "line", x0 = 1, x1 = 1, y0 = 0, y1 = 1,
                           yref = "paper",
                           line = list(color = "red", dash = "dash", width = 1))),
        legend = list(title = list(text = "Cancer")),
        margin = list(l = 160),
        paper_bgcolor = bl$paper_bgcolor,
        plot_bgcolor  = bl$plot_bgcolor,
        font          = bl$font
      )
  })

  # Colocalization table
  output$coloc_table <- renderDT({
    req(master_ev)
    master_ev %>%
      select(protein, cancer_mr, tier_short,
             coloc_PPH4_abf, coloc_PPH4_susie, coloc_verdict) %>%
      arrange(tier_short, protein) %>%
      datatable(
        rownames = FALSE, class = "cell-border stripe hover",
        colnames = c("Protein", "Cancer", "Tier",
                     "PPH4 (ABF)", "PPH4 (SuSiE)", "Verdict"),
        options  = list(pageLength = 20, dom = "frtip", autoWidth = TRUE)
      ) %>%
      formatRound(c("coloc_PPH4_abf", "coloc_PPH4_susie"), digits = 3) %>%
      formatStyle("tier_short",
                  backgroundColor = styleEqual(names(tier_pal), unname(tier_pal)))
  })

  # ER concordance scatter
  output$er_scatter <- renderPlotly({
    req(master_ev)
    er  <- master_ev %>% filter(!is.na(or_ERpos) & !is.na(or_ERneg))
    rng <- range(c(er$or_ERpos, er$or_ERneg), na.rm = TRUE)
    pad <- diff(rng) * 0.08
    lo  <- rng[1] - pad; hi <- rng[2] + pad

    bl <- pl_layout()
    plot_ly(er, x = ~or_ERpos, y = ~or_ERneg,
            type = "scatter", mode = "markers+text",
            text = ~protein, textposition = "top center",
            hovertext = ~paste0("<b>", protein, "</b><br>",
                                "ER+: ", round(or_ERpos, 3),
                                " (p=", signif(pval_ERpos, 2), ")<br>",
                                "ER−: ", round(or_ERneg, 3),
                                " (p=", signif(pval_ERneg, 2), ")"),
            hoverinfo = "text",
            marker = list(size = 11, color = "#005b96", opacity = 0.85)) %>%
      layout(
        xaxis = list(title = "OR — ER positive", range = c(lo, hi),
                     gridcolor = bl$gridcolor),
        yaxis = list(title = "OR — ER negative", range = c(lo, hi),
                     gridcolor = bl$gridcolor),
        shapes = list(
          list(type = "line", x0 = lo, x1 = hi, y0 = lo, y1 = hi,
               line = list(color = "grey", dash = "dot")),
          list(type = "line", x0 = 1, x1 = 1, y0 = lo, y1 = hi,
               line = list(color = "#aaa", dash = "dash", width = 1)),
          list(type = "line", x0 = lo, x1 = hi, y0 = 1, y1 = 1,
               line = list(color = "#aaa", dash = "dash", width = 1))
        ),
        paper_bgcolor = bl$paper_bgcolor,
        plot_bgcolor  = bl$plot_bgcolor,
        font          = bl$font
      )
  })

  # Evidence table
  output$evidence_table <- renderDT({
    req(master_ev)
    master_ev %>%
      select(protein, cancer_mr, tier_short,
             mr_or, mr_or_lo, mr_or_hi, mr_pval, mr_fdr,
             coloc_PPH4_best, coloc_verdict, magma_breast_p) %>%
      arrange(tier_short, cancer_mr) %>%
      datatable(
        extensions = "Buttons", rownames = FALSE,
        class = "cell-border stripe hover",
        colnames = c("Protein", "Cancer", "Tier",
                     "OR", "OR low", "OR high", "P-value", "FDR",
                     "PPH4 (best)", "Coloc verdict", "MAGMA p"),
        options = list(
          pageLength = 20, dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          autoWidth = TRUE, scrollX = TRUE
        )
      ) %>%
      formatRound(c("mr_or", "mr_or_lo", "mr_or_hi", "coloc_PPH4_best"), digits = 3) %>%
      formatSignif(c("mr_pval", "mr_fdr", "magma_breast_p"), digits = 3) %>%
      formatStyle("tier_short",
                  backgroundColor = styleEqual(names(tier_pal), unname(tier_pal)))
  })

  # Druggability table
  output$druggability_table <- renderDT({
    req(druggability)
    druggability %>%
      select(protein, approved_name, n_known_drugs,
             tractability_SM, tractability_AB,
             top_drug, top_drug_phase, top_indication) %>%
      datatable(
        extensions = "Buttons", rownames = FALSE,
        class = "cell-border stripe hover",
        colnames = c("Protein", "Gene name", "Known drugs",
                     "SM tractability", "AB tractability",
                     "Top drug", "Phase", "Indication"),
        options = list(
          pageLength = 20, dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          autoWidth = TRUE, scrollX = TRUE
        )
      )
  })

  # Mediation table
  output$mediation_table <- renderDT({
    req(mediation)
    mediation %>%
      select(protein, metabolite, cancer,
             p_indirect, prop_med_pct, coloc_evidence_class) %>%
      datatable(
        extensions = "Buttons", rownames = FALSE,
        class = "cell-border stripe hover",
        colnames = c("Protein", "Metabolite", "Cancer",
                     "Indirect p", "% Mediated", "Coloc class"),
        options = list(
          pageLength = 20, dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          autoWidth = TRUE, scrollX = TRUE
        )
      ) %>%
      formatSignif("p_indirect", digits = 3) %>%
      formatRound("prop_med_pct", digits = 1)
  })
}

shinyApp(ui, server)
