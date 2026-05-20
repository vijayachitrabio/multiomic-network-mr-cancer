#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
})

out_fig <- "results/figures"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

box <- function(label, x, y, w, h, gp = gpar(fill = "white", col = "#1F1F1F", lwd = 2),
                text_gp = gpar(col = "#1F1F1F", fontsize = 10.5), r = 0.018) {
  grid.roundrect(x, y, w, h, r = unit(r, "snpc"), gp = gp)
  grid.text(label, x, y, gp = text_gp)
}

draw_arrow <- function(x0, y0, x1, y1, lwd = 1.8) {
  grid.segments(
    x0, y0, x1, y1,
    gp = gpar(col = "#1F1F1F", lwd = lwd),
    arrow = grid::arrow(length = unit(0.12, "in"), type = "closed")
  )
}

draw_panel <- function() {
  grid.newpage()

  # Panel headings
  grid.text("Trial analogy", x = 0.25, y = 0.910, gp = gpar(fontsize = 12.5, fontface = "bold"))
  grid.text("This study: gene-based protein analysis", x = 0.74, y = 0.910, gp = gpar(fontsize = 12.5, fontface = "bold"))

  # Left panel: RCT analogy
  box("randomization\nby investigator", 0.25, 0.835, 0.38, 0.062)
  box("control", 0.145, 0.725, 0.16, 0.085)
  box("treatment", 0.355, 0.725, 0.16, 0.085)
  box("exposure\nlower", 0.145, 0.580, 0.18, 0.078)
  box("exposure\nhigher", 0.355, 0.580, 0.18, 0.078)
  box("outcome\ncompared", 0.25, 0.410, 0.38, 0.075)

  draw_arrow(0.25, 0.803, 0.145, 0.770)
  draw_arrow(0.25, 0.803, 0.355, 0.770)
  draw_arrow(0.145, 0.682, 0.145, 0.620)
  draw_arrow(0.355, 0.682, 0.355, 0.620)
  draw_arrow(0.145, 0.541, 0.205, 0.450)
  draw_arrow(0.355, 0.541, 0.295, 0.450)

  box(
    "measured and unmeasured\nrisk factors balanced by design",
    0.25, 0.500, 0.34, 0.060,
    gp = gpar(fill = "#F7F8F9", col = "#1F1F1F", lwd = 1.6),
    text_gp = gpar(col = "#1F1F1F", fontsize = 8.5)
  )

  # Right panel: our MR framework
  box("individuals grouped\nby inherited gene variant", 0.74, 0.835, 0.38, 0.062)
  box("protein-changing\nvariant absent", 0.625, 0.725, 0.17, 0.092)
  box("protein-changing\nvariant present", 0.855, 0.725, 0.17, 0.092)
  box("protein level\nlower", 0.625, 0.580, 0.18, 0.078)
  box("protein level\nhigher", 0.855, 0.580, 0.18, 0.078)

  box("possible\nmetabolite pathway", 0.740, 0.450, 0.26, 0.078)
  box("breast or\nendometrial cancer risk", 0.740, 0.290, 0.35, 0.078)

  draw_arrow(0.74, 0.803, 0.625, 0.773)
  draw_arrow(0.74, 0.803, 0.855, 0.773)
  draw_arrow(0.625, 0.679, 0.625, 0.620)
  draw_arrow(0.855, 0.679, 0.855, 0.620)
  draw_arrow(0.625, 0.541, 0.695, 0.490)
  draw_arrow(0.855, 0.541, 0.785, 0.490)
  draw_arrow(0.740, 0.410, 0.740, 0.333)
  draw_arrow(0.625, 0.541, 0.680, 0.333)
  draw_arrow(0.855, 0.541, 0.800, 0.333)

  box(
    "checks used:\ninstrument strength, direction tests,\nshared-signal tests, replication, sensitivity analyses",
    0.740, 0.180, 0.43, 0.085,
    gp = gpar(fill = "#F7F8F9", col = "#1F1F1F", lwd = 1.6),
    text_gp = gpar(col = "#1F1F1F", fontsize = 8.4)
  )
  draw_arrow(0.740, 0.251, 0.740, 0.223, lwd = 1.4)

  # Divider
  grid.segments(0.50, 0.14, 0.50, 0.90, gp = gpar(col = "#D3D6D9", lwd = 1.2))
}

png(file.path(out_fig, "sfig_mr_design_multiomic.png"),
    width = 10.8, height = 7.0, units = "in", res = 340, bg = "white")
draw_panel()
dev.off()

pdf(file.path(out_fig, "sfig_mr_design_multiomic.pdf"),
    width = 10.8, height = 7.0, bg = "white")
draw_panel()
dev.off()

cat("Saved supplementary MR design figure:\n")
cat(" - results/figures/sfig_mr_design_multiomic.png\n")
cat(" - results/figures/sfig_mr_design_multiomic.pdf\n")
