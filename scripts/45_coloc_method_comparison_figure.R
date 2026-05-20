#!/usr/bin/env Rscript
## Script 45 (v2): coloc.abf vs coloc.susie comparison figure — clean, no overlaps
## Output: results/figures/fig_coloc_method_comparison.png/.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})

proj    <- "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/MultiOmic_Network_MR_Project"
dat_f   <- file.path(proj, "results/tables/STable8_protein_coloc.csv")
out_dir <- file.path(proj, "results/figures")

d <- fread(dat_f)
d[is.na(PPH4_abf),   PPH4_abf   := 0]
d[is.na(PPH4_susie), PPH4_susie := 0]

# Quadrant classification (short labels for legend)
d[, quadrant := fcase(
  PPH4_susie >= 0.8 & PPH4_abf < 0.5,  "SuSiE resolves (ABF missed)",
  PPH4_susie >= 0.8 & PPH4_abf >= 0.5, "Both agree: strong",
  PPH4_abf   >= 0.5 & PPH4_susie < 0.5,"ABF only",
  default = "Neither"
)]

col_map <- c(
  "SuSiE resolves (ABF missed)" = "#C0392B",
  "Both agree: strong"          = "#1A5276",
  "ABF only"                    = "#E67E22",
  "Neither"                     = "#AAB7B8"
)
shape_map <- c(
  "SuSiE resolves (ABF missed)" = 17,   # triangle
  "Both agree: strong"          = 16,   # circle
  "ABF only"                    = 15,   # square
  "Neither"                     = 16
)

# No manual nudges — let ggrepel spread naturally within expanded canvas

# Zone labels: anchored to the shaded-region corners, far from data points
zone_df <- data.frame(
  x     = c(0.24,  0.38,  0.82,  0.82),
  y     = c(0.83,  0.05,  0.83,  0.05),
  hjust = c(0.5,   0.5,   0,     0),
  vjust = c(0,     0,     0,     0),
  col   = c("#C0392B","#9B9B9B","#1A5276","#E67E22"),
  label = c("SuSiE resolves\n(ABF missed)", "Neither", "Both\nstrong", "ABF\nonly")
)

p <- ggplot(d, aes(x = PPH4_abf, y = PPH4_susie,
                   colour = quadrant, shape = quadrant)) +
  # shaded quadrants
  annotate("rect", xmin=-0.02, xmax=0.5,  ymin=0.8,  ymax=1.04,
           fill="#C0392B", alpha=0.07) +
  annotate("rect", xmin=0.8,  xmax=1.04, ymin=0.8,  ymax=1.04,
           fill="#1A5276", alpha=0.07) +
  # threshold + diagonal
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey65", linewidth=0.4) +
  geom_hline(yintercept=0.8, linetype="dotted", colour="grey40", linewidth=0.35) +
  geom_vline(xintercept=0.8, linetype="dotted", colour="grey40", linewidth=0.35) +
  # points
  geom_point(size=6, alpha=0.93) +
  # labels — ggrepel handles spreading; force_pull=0 prevents collapse back to point
  geom_text_repel(
    aes(label = protein),
    fontface           = "italic",
    size               = 4.0,
    colour             = "grey10",
    segment.color      = "grey55",
    segment.size       = 0.35,
    min.segment.length = 0.05,
    box.padding        = 0.55,
    point.padding      = 0.4,
    force              = 6,
    force_pull         = 0.5,
    max.overlaps       = 20,
    xlim               = c(-0.12, 1.12),
    ylim               = c(-0.12, 1.12),
    seed               = 99,
    show.legend        = FALSE
  ) +
  # tiny quadrant corner labels
  geom_text(data=zone_df,
            aes(x=x, y=y, label=label, hjust=hjust, vjust=vjust, colour=NULL, shape=NULL),
            colour=zone_df$col, size=3.0, fontface="bold", lineheight=0.85,
            inherit.aes=FALSE) +
  scale_colour_manual(values=col_map, name=NULL) +
  scale_shape_manual(values=shape_map,  name=NULL) +
  scale_x_continuous(limits=c(-0.12, 1.12), breaks=seq(0,1,0.25)) +
  scale_y_continuous(limits=c(-0.12, 1.12), breaks=seq(0,1,0.25)) +
  labs(
    x       = "coloc.abf  PPH4",
    y       = "coloc.susie  PPH4  (SuSiE fine-mapping + LD)",
    title   = "coloc.abf  vs  coloc.susie — colocalization posteriors",
    subtitle= "8 priority proteins  ·  breast cancer GWAS (N = 228,951)  ·  FinnGen pQTL (N = 619)",
    caption = "Dotted lines: PPH4 = 0.80 threshold.  Dashed diagonal: perfect agreement between methods."
  ) +
  theme_bw(base_size=13) +
  theme(
    plot.title      = element_text(face="bold", size=13.5, colour="#1A5276"),
    plot.subtitle   = element_text(size=10, colour="grey35"),
    plot.caption    = element_text(size=8.5, colour="grey55"),
    legend.position = "bottom",
    legend.key.size = unit(0.45,"cm"),
    legend.text     = element_text(size=10),
    panel.grid.minor= element_blank(),
    panel.grid.major= element_line(colour="grey92"),
    axis.title      = element_text(size=11.5),
    plot.margin     = margin(10, 14, 8, 10)
  ) +
  guides(colour=guide_legend(override.aes=list(size=4)),
         shape =guide_legend(override.aes=list(size=4)))

ggsave(file.path(out_dir, "fig_coloc_method_comparison.png"),
       p, width=8, height=7.2, dpi=300, bg="white")
ggsave(file.path(out_dir, "fig_coloc_method_comparison.pdf"),
       p, width=8, height=7.2)

message("✓ fig_coloc_method_comparison.png/.pdf saved")
